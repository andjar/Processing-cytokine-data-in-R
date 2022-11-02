InpregObject <- R6::R6Class("InpregObject",
      lock_objects = FALSE,
      public = list(
        # Set up
        filenames = list(
          FI = NA,
          biorad = NA,
          expected = NA
        ),
        get_all = FALSE,
        use_FI = TRUE,
        use_log = TRUE,
        reference_batch = NULL,

        below_OOR = -1,
        above_OOR = -2,

        norm_function = NULL,

        # Logging
        log_file = NULL,
        logger = NULL,
        log_level = NULL,
        do_debug = FALSE,

        # Data
        dfs = list(
          FI = data.table(),
          biorad = data.table(),
          expected = data.table()
        ),
        dfs_missing = list(
          FI = data.table(),
          biorad = data.table()
        ),
        res.drm = NULL,

        # Plotting
        my_theme = ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom"),

        # Functions
        initialize = function(...) {

          # Set properties
          inputs <- list(...)
          for (i in seq_along(inputs)) {
            self[[names(inputs)[i]]] <- inputs[[i]]
          }

          # Make log
          self$log_file <- paste0(format(Sys.time(), format = "%Y%m%dT%H%M%S"), ".log")
          self$log_level <- ifelse(self$do_debug, "DEBUG", "INFO")
          self$logger <- log4r::logger(self$log_level, appenders = list(log4r::console_appender(), log4r::file_appender(self$log_file)))

          # Import data
          self$dfs[["FI"]] <- import_data(filename = self$filenames[["FI"]], get_all = self$get_all)
          self$dfs[["biorad"]] <- import_data(filename = self$filenames[["biorad"]], get_all = self$get_all)
          self$dfs[["expected"]] <- fread(self$filenames[["expected"]])
          self$dfs[["expected"]][, batch := factor(batch)]
          self$dfs[["expected"]] <- melt(self$dfs[["expected"]], id.vars = c("study", "batch", "type"), value.name = "concentration")
          self$log("Files imported")


          self$dfs_missing[["FI"]] <- self$wide_to_long(self$df(FI = TRUE))[value < 0]
          self$dfs_missing[["biorad"]] <- self$wide_to_long(self$df(FI = FALSE))[value < 0]

          for (i in self$cytokines) {
            set(self$dfs[["FI"]], i = self$dfs[["FI"]][get(i) < 0, which = TRUE], j = i, value = NA)
            set(self$dfs[["biorad"]], i = self$dfs[["biorad"]][get(i) < 0, which = TRUE], j = i, value = NA)
          }

          if (self$use_FI) {
            self$log("Using FI")
          }

          if (self$use_log) {
            self$log("Using log_e-values")
            for (i in self$cytokines) {
              set(self$dfs[["FI"]], i = self$dfs[["FI"]][get(i) >= 0, which = TRUE], j = i, value = log(self$dfs[["FI"]][get(i) >= 0, get(i)]))
              set(self$dfs[["biorad"]], i = self$dfs[["biorad"]][get(i) >= 0, which = TRUE], j = i, value = log(self$dfs[["biorad"]][get(i) >= 0, get(i)]))
            }
          }

          private$dfs0 <- copy(self$dfs)
          private$dfs_missing0 <- copy(self$dfs_missing)

          # https://stackoverflow.com/questions/39914775/updating-method-definitions-in-r6-object-instance
          self$norm_function <- function(x) mean(x, na.rm = TRUE)

          },
        df = function(FI = NULL, types = NULL, batches = NULL) {
          if(is.null(FI)) {
            FI <- self$use_FI
          }
          if (FI) {
            dfo <- self$dfs[["FI"]]
          } else {
            dfo <- self$dfs[["biorad"]]
          }

          if (!is.null(types)) {
            dfo <- dfo[substr(type, 1, 1) %in% types]
          }

          if (!is.null(batches)) {
            dfo <- dfo[batch %in% batches]
          } else {
            if (!is.null(private$batches)) {
              dfo <- dfo[batch %in% private$batches]
            }
          }

          if (!is.null(private$exclude_variables)) {
            dfo <- dfo[, .SD, .SDcols = colnames(dfo)[!colnames(dfo) %in% private$exclude_variables]]
          }

          return(dfo)
        },
        df_missing = function() {
          self$log("Showing missing on biorad data")
          dfo <- self$dfs_missing[["biorad"]]
          if (!is.null(private$types)) {
            dfo <- dfo[substr(type, 1, 1) %in% private$types]
          }
          if (!is.null(private$batches)) {
            dfo <- dfo[batch %in% private$batches]
          }
          if (!is.null(private$exclude_variables)) {
            dfo <- dfo[, .SD, .SDcols = colnames(dfo)[!colnames(dfo) %in% private$exclude_variables]]
          }
          return(dfo)
        },
        wide_to_long = function(df = NULL) {
          if (is.null(df)) {
            df <- self$df()
          }

          df <- melt(df, id.vars = colnames(df)[!colnames(df) %in% self$cytokines])
          df[, value := as.numeric(value)]

          return(df)
        },
        pca = function(types = NULL) {
          df <- self$df(types = types)
          df_plot <- df[, .SD, .SDcols = self$cytokines]
          df <- df[, .SD, .SDcols = colnames(df)[!colnames(df) %in% colnames(df_plot)]]
          df <- df[complete.cases(df_plot), ]
          n1 <- nrow(df_plot)
          df_plot <- na.omit(df_plot)
          n2 <- nrow(df_plot)
          if (n2 != n1) {
            self$log(paste0("Removed ", n2-n1, " rows due to missing values"), level = "WARN")
          }
          df_pca  <- prcomp(df_plot, scale = TRUE, center = TRUE)
          df_score <- as.data.table(df_pca$x)
          df_score <- cbind(df, df_score)
          return(df_score)
        },
        log = function(message, level = "INFO") {
          log4r::levellog(self$logger, level = level, message)
        },
        set_type = function(type = NULL) {
          private$types <- type
          self$log(paste0("Only using samples of type ", type), level = "WARN")
        },
        set_batch = function(batch = NULL) {
          private$batches <- batch
          self$log(paste0("Only using samples from batch(es) ", batch), level = "WARN")
        },
        rm_variables = function(variables = NULL) {
          private$exclude_variables <- variables
          self$log(paste0("Excluding ", variables), level = "WARN")
        },
        plot = function(type = "pca") {
          self$plot_pca()
        },
        plot_pca = function(components = c(1, 2), color = "batch", highlight = NULL) {
          df <- self$pca(types = private$types)
          g <- ggplot(df,
                      aes_string(
                        x = paste0("PC", components[1]),
                        y = paste0("PC", components[2])
                      )
          )

          if (is.null(color) & is.null(highlight)) {
            g <- g + geom_point()
          } else {
            if (is.null(highlight)) {
              g <- g + geom_point(aes_string(color = color), alpha = 0.5) + labs(color = color)
            } else {
              g <- g + geom_point(aes(color = ifelse(bp %in% highlight, bp, "Other")), alpha = 0.5) + labs(color = "Selected plates")
            }
          }

          g <- g + self$my_theme
          return(g)
        },
        plot_boxplot = function(x_axis = "batch") {
          df <- self$wide_to_long(self$df(types = private$types))
          ggplot(df, aes_string(x_axis, "value")) +
            geom_boxplot() +
            facet_wrap("variable", scales = "free_y") +
            self$my_theme
        },
        plot_missing = function(x_axis = "bp", prop = TRUE) {
          df <- self$df_missing()
          df_plot <- df[, .(too_low = sum(value == self$below_OOR),
                       too_high = sum(value == self$above_OOR)), by = c("variable", x_axis)]
          df_plot[, cytokine := variable]
          df_plot[, variable := NULL]
          df_plot <- melt(df_plot, id.vars = c(x_axis, "cytokine"))
          if (prop) {
            self$log("Showing proportion")
            df_all <- self$df()[, .N, by = x_axis]
            df_plot <- merge(df_plot, df_all)
            df_plot[, value := value / N]
          }
          ggplot(df_plot, aes_string(x_axis, "value")) +
            geom_point(aes_string(shape = "variable")) +
            facet_wrap("cytokine") +
            coord_flip() +
            self$my_theme
        },
        plot_qc = function(types = c("blank", "internkontroll", "serumkontroll", "internstandard"), scales = "fixed") {
          df <- self$wide_to_long(self$df(types = "C"))
          df <- df[studyno %in% types, ]
          ggplot(df, aes(bp, value, shape = studyno)) +
            geom_point() +
            facet_wrap(~variable, scales = scales) +
            coord_flip() +
            self$my_theme
        },
        plot_time = function(scales = "fixed") {
          df <- self$wide_to_long(self$df(types = "X"))

          ggplot(df, aes(time, value, group = studyno)) +
            geom_point() +
            geom_line() +
            facet_wrap(~variable, scales = scales) +
            self$my_theme
        },
        plot_medians = function(x_axis = "bp", types = NULL, scales = "fixed") {
          df <- self$wide_to_long(self$df(types = private$types))
          df[, sample_type := substr(type, 1, 1)]
          df[, sample_type := ifelse(sample_type == "C", studyno, sample_type)]
          if (!is.null(types)) {
            df <- df[sample_type %in% types, ]
          }
          df <- df[, .(m = median(value, na.rm = TRUE), l = quantile(value, 1/4, na.rm = TRUE), h = quantile(value, 3/4, na.rm = TRUE)), by = c("variable", "sample_type", x_axis)]
          ggplot(df, aes_string(x = x_axis, y = "m", ymin = "l", ymax = "h", shape = "sample_type")) +
            geom_pointrange(position = position_dodge2(width = 0.5)) +
            facet_wrap(~variable, scales = scales) +
            self$my_theme + coord_flip()
        },
        plot_standards = function(scales = "free_x", highlight = NULL) {
          df <- self$wide_to_long(self$df(types = "S"))
          df <- merge(df, self$dfs[["expected"]])
          if (is.null(highlight)) {
            ggplot(df, aes(concentration, value, group = bp)) +
              geom_point() + geom_line() +
              facet_wrap(~variable, scales = scales) +
              scale_x_log10() +
              self$my_theme
          } else {
            ggplot(df, aes(concentration, value, group = bp, color = ifelse(bp %in% highlight, bp, "Other"))) +
              geom_point() + geom_line() +
              facet_wrap(~variable, scales = scales) +
              scale_x_log10() +
              self$my_theme + labs(color = "Selected")
          }

        },
        find_duplicates = function(types = "X") {
          df <- self$df(types = ifelse(is.null(types), private$types, types))
          df[, study_id_time := paste(study, studyno, time)]
          dup_ids <- df[duplicated(study_id_time), unique(study_id_time)]
          df2 <- df[study_id_time %in% dup_ids, ]
          return(df2)
        },
        find_duplicates_2 = function(ref = "2", types = "X") {
          df <- self$find_duplicates(types = types)
          df_ref <- df[batch == ref, ]
          df <- df[batch != ref, ]
          df <- self$wide_to_long(df)
          df_ref <- self$wide_to_long(df_ref)
          df <- merge(df_ref[, .SD, .SDcols = c("study_id_time", "variable", "batch", "bp", "value")],
                      df[, .SD, .SDcols = c("study_id_time", "variable", "batch", "bp", "value")],
                      by = c("study_id_time", "variable"),
                      all.y = FALSE,
                      all.x = TRUE)
        },
        get_duplicate_factor = function(ref = "2", types = "X") {
          df <- self$find_duplicates_2(ref = ref, types = types)
          # x is for the reference batch, y for the others
          return(
            df[, .(mref = self$norm_function(value.x), my = self$norm_function(value.y)), by = c("variable", "batch.y")][, f := mref / my]
          )
        },
        plot_duplicates = function(ref = "2", scales = "free") {
          df <- self$find_duplicates_2(ref = ref)
          ggplot(df, aes(value.x, value.y, color = bp.y)) +
            geom_point() + geom_smooth(se = FALSE, method = "lm") +
            facet_wrap(~variable, scales = scales) +
            labs(x = paste("Batch", ref)) +
            self$my_theme
        },
        get_plate_factor = function(types = NULL, by = NULL, ref_study = NULL, ref_time = NULL, ref_plates = NULL) {
          if (is.null(by)) {
            if (is.null(ref_study)) {
              if (is.null(ref_time)) {
                by = c("batch", "plate")
              } else {
                by = c("batch", "plate", "time")
              }
            } else {
              if (is.null(ref_time)) {
                by = c("batch", "plate", "study")
              } else {
                by = c("batch", "plate", "study", "time")
              }
            }
          }

          df <- self$wide_to_long(
            self$df(types = ifelse(is.null(types), private$types, types))
            )[, .(m = self$norm_function(value)), by = c(by, "variable")]

          if (is.null(ref_study)) {

            if (is.null(ref_time)) {

              if (is.null(ref_plates)) {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[, .(m = self$norm_function(value)), by = c("batch", "variable")]

              } else {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[actualPlate %in% ref_plates, .(m = self$norm_function(value)), by = c("batch", "variable")]

              }

              # x refers to grand mean, y to local mean
              df <- merge(df_ref, df, by = c("batch", "variable"))
              df[, f := m.x / m.y]

            } else {

              if (is.null(ref_plates)) {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[, .(m = self$norm_function(value)), by = c("batch", "time", "variable")]

              } else {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[actualPlate %in% ref_plates, .(m = self$norm_function(value)), by = c("batch", "time", "variable")]

              }

              df_ref <- df_ref[study %in% ref_time, ]
              df <- df[study %in% ref_time, ]

              # x refers to grand mean, y to local mean
              df <- merge(df_ref, df, by = c("batch", "time", "variable"))
              df[, f := m.x / m.y]

            }

          } else {

            if (is.null(ref_time)) {

              if (is.null(ref_plates)) {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[, .(m = self$norm_function(value)), by = c("batch", "study", "variable")]
                df_ref <- df_ref[study == ref_study, ]
                df <- df[study == ref_study, ]

              } else {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[actualPlate %in% ref_plates, .(m = self$norm_function(value)), by = c("batch", "study", "variable")]
                df_ref <- df_ref[study == ref_study, ]
                df <- df[study == ref_study, ]

              }

              # x refers to grand mean, y to local mean
              df <- merge(df_ref, df, by = c("batch", "study", "variable"))
              df[, f := m.x / m.y]

            } else {



              if (is.null(ref_plates)) {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[, .(m = self$norm_function(value)), by = c("batch", "time", "study", "variable")]
                df_ref <- df_ref[study == ref_study & time %in% ref_time, ]
                df <- df[study == ref_study & time %in% ref_time, ]

              } else {

                df_ref <- self$wide_to_long(
                  self$df(types = ifelse(is.null(types), private$types, types))
                )[actualPlate %in% ref_plates, .(m = self$norm_function(value)), by = c("batch", "time", "study", "variable")]
                df_ref <- df_ref[study == ref_study & time %in% ref_time, ]
                df <- df[study == ref_study & time %in% ref_time, ]

              }

              # x refers to grand mean, y to local mean
              df <- merge(df_ref, df, by = c("batch", "time", "study", "variable"))
              df[, f := m.x / m.y]

            }

          }


          return(df)
        },
        adjust_plate_effect = function(df = NULL, types = NULL, ref_study = NULL, ref_time = NULL) {
          self$log("Adjusting both FI and biorad values")
          if (is.null(df)) {
            df <- self$get_plate_factor(types = types, ref_study = ref_study, ref_time = ref_time)
          }
          for (b in unique(df$batch)) {
            for (p in unique(df$plate[df$batch == b])) {
              for (v in self$cytokines) {
                set(self$dfs[["FI"]],
                    i = self$dfs[["FI"]][batch == b & plate == p, which = TRUE],
                    j = v,
                    value = self$dfs[["FI"]][batch == b & plate == p, get(v)] * df[batch == b & plate == p & variable == v, f]
                    )
                set(self$dfs[["biorad"]],
                    i = self$dfs[["biorad"]][batch == b & plate == p, which = TRUE],
                    j = v,
                    value = self$dfs[["biorad"]][batch == b & plate == p, get(v)] * df[batch == b & plate == p & variable == v, f]
                )
              }
            }
          }
        },
        adjust_plate_effect_by_time = function(df = NULL) {
          self$log("Adjusting both FI and biorad values")

          for (b in unique(df$batch)) {
            for (p in unique(df$plate[df$batch == b])) {
              for (t in unique(df$time[df$batch == b & df$plate == p])) {
                for (v in self$cytokines) {
                  set(self$dfs[["FI"]],
                      i = self$dfs[["FI"]][batch == b & plate == p & time == t, which = TRUE],
                      j = v,
                      value = self$dfs[["FI"]][batch == b & plate == p & time == t, get(v)] * df[batch == b & plate == p & time == t & variable == v, f]
                  )
                  set(self$dfs[["biorad"]],
                      i = self$dfs[["biorad"]][batch == b & plate == p & time == t, which = TRUE],
                      j = v,
                      value = self$dfs[["biorad"]][batch == b & plate == p & time == t, get(v)] * df[batch == b & plate == p & time == t & variable == v, f]
                  )
                }
              }
            }
          }
        },
        adjust_batch_effect = function(ref = "2", types = "X") {

          self$log("Adjusting both FI and biorad values")
          df <- self$get_duplicate_factor(ref = ref, types = types)
          for (b in unique(df$batch.y)) {
            if (b != ref) {
              for (v in self$cytokines) {
                set(self$dfs[["FI"]],
                    i = self$dfs[["FI"]][batch == b, which = TRUE],
                    j = v,
                    value = self$dfs[["FI"]][batch == b, get(v)] * df[batch.y == b & variable == v, f]
                )
                set(self$dfs[["biorad"]],
                    i = self$dfs[["biorad"]][batch == b, which = TRUE],
                    j = v,
                    value = self$dfs[["biorad"]][batch == b, get(v)] * df[batch.y == b & variable == v, f]
                )
              }
            }
          }
        },
        df_raw = function(FI = NULL, types = NULL, batches = NULL) {
          if(is.null(FI)) {
            FI <- self$use_FI
          }
          if (FI) {
            dfo <- private$dfs0[["FI"]]
          } else {
            dfo <- private$dfs0[["biorad"]]
          }

          if (!is.null(types)) {
            dfo <- dfo[substr(type, 1, 1) %in% types]
          }

          if (!is.null(batches)) {
            dfo <- dfo[batch %in% batches]
          } else {
            if (!is.null(private$batches)) {
              dfo <- dfo[batch %in% private$batches]
            }
          }

          if (!is.null(private$exclude_variables)) {
            dfo <- dfo[, .SD, .SDcols = colnames(dfo)[!colnames(dfo) %in% private$exclude_variables]]
          }

          return(dfo)
        },
        df_cmp = function(FI = NULL, types = NULL, batches = NULL, wide = TRUE) {

          if (wide) {
            df <- self$wide_to_long(self$df(FI = FI, types = types, batches = batches))
            df_r <- self$wide_to_long(self$df_raw(FI = FI, types = types, batches = batches))

            colnames(df)[colnames(df) == "value"] <- "value.adj"
            colnames(df_r)[colnames(df_r) == "value"] <- "value.raw"
            merge_by <- colnames(df)
            merge_by <- merge_by[merge_by != "value.adj"]

            return(
              merge(df, df_r, by = merge_by)
            )
          } else {
            df <- self$wide_to_long(self$df(FI = FI, types = types, batches = batches))
            df_r <- self$wide_to_long(self$df_raw(FI = FI, types = types, batches = batches))
            df[ , adjusted := TRUE]
            df_r[ , adjusted := FALSE]

            return(
              rbind(df, df_r)
            )
          }

        },
        c_o_v = function(df, by = c("batch", "plate")) {
          return(
            df[, .(cov = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)), by = c(by, "variable")]
          )
        },
        get_concentrations = function(ref = "2") {
          df <- self$wide_to_long(self$df(FI = TRUE))
          df <- merge(df, self$dfs[["expected"]][batch == ref,], by = c("batch", "study", "type", "variable"), all.x = TRUE, all.y = FALSE)
          if (!self$use_log) {
            df[, value := log(value)]
          }

          df[, well_role := ifelse(substr(type, 1, 1) == "S" & batch == ref, "Standard", "Unknown")]
          df[, analyte := variable]
          df[, assay_id := ""]
          df[, dilution := ifelse(substr(type, 1, 1) == "X", 4, 1)]
          df[, sample_id := ifelse(substr(type, 1, 1) == "S", type, unique_sample_ID)]

          self$log("Starting nCal - this will take some time...")
          self$res.drm <- nCal::ncal(value~concentration, df,
                        return.fits = TRUE,
                        unk.replicate = 1,
                        bcrm.fit = FALSE,
                        force.fit = TRUE,
                        var.model = "power",
                        find.LOD = TRUE,
                        find.LOQ = TRUE,
                        find.best.dilution = FALSE,
                        plot = FALSE,
                        plot.se.profile=TRUE)
          self$log("--> Finished")

          self$dfs[["LOD"]] <- self$get_attr(attr = "LOD")
          self$dfs[["LOQ"]] <- self$get_attr(attr = "LOQ")
          self$res.drm <- as.data.table(self$res.drm)
          self$res.drm[, value := ifelse(self$use_log, est.log.conc, est.conc)]
          self$dfs[["FIr"]] <- self$dfs[["FI"]]
          self$dfs[["FI"]] <- dcast(data = as.data.table(res.drm),
                                    formula(paste0(paste(colnames(res.drm)[colnames(res.drm) %in% colnames(self$dfs[["FI"]])], collapse = "+"), "~variable")),
                                    value.var = "value")

         },
        get_attr = function(attr = NULL) {
          return(attr(self$res.drm, attr))
        }
      ),
      active = list(
        cytokines = function() {
          cytokines <- c("IL_1b", "IL_1ra", "IL_2", "IL_4", "IL_5", "IL_6", "IL_7", "IL_8", "IL_9", "IL_10", "IL_12", "IL_13", "IL_15", "IL_17", "eotaxin", "FGF_b", "G_CSF", "GM_CSF", "IFN_g", "IP_10", "MCP_1", "MIP_1a", "PDGF_bb", "MIP_1b", "RANTES", "TNF_a", "VEGF")
          cytokines <- cytokines[!cytokines %in% private$exclude_variables]
        }
      ),
      private = list(
        batches = NULL,
        types = NULL,
        exclude_variables = NULL,
        dfs0 = NULL,
        dfs_missing0 = NULL
      )
)
