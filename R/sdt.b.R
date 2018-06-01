
sdtClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "sdtClass",
    inherit = sdtBase,
    private = list(
        #### Init + run functions ----
        .init = function() {

            private$.initSdtGroupTable()
            private$.initSdtSubjTable()

        },
        .run = function() {

            if (is.null(self$options$res) || is.null(self$options$stim))
                return()

            data <- private$.cleanData()
            results <- private$.compute(data)

            sumStats <- private$.sumStats(data, results)

            private$.populateSdtGroupTable(sumStats$group)
            private$.populateSdtSubjTable(sumStats$subjs)

            private$.preparedPrimePlot(sumStats)
        },

        #### Compute results ----
        .compute = function(data) {

            path <- system.file("models/sdt.txt", package = "resWagner")

            # Sampling Parameters
            nChains = self$options$nChains
            nBurnin = self$options$nBurnin
            nSamples = self$options$nSamples
            nThin = self$options$nThin

            # Parameters to be monitored
            parameters <- c('dMu', 'cMu', 'd', 'c', 'dTypSubj', 'cTypSubj')

            # Draw posterior samples
            model <- rjags::jags.model(file=path, data=data, n.chains=nChains)

            if (nBurnin > 0)
                update(model, nBurnin)

            samples <- rjags::jags.samples(model, parameters, thin=nThin, n.iter=nSamples)

            r <- list(samples=samples)

            return(r)
        },

        #### Init tables/plots functions ----
        .initSdtGroupTable = function() {

            table <- self$results$sdtGroup

            ciWidth <- self$options$ciWidth

            table$getColumn('dPrimeLower')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))
            table$getColumn('dPrimeUpper')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))
            table$getColumn('cLower')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))
            table$getColumn('cUpper')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))

        },
        .initSdtSubjTable = function() {

            table <- self$results$sdtSubj

            ciWidth <- self$options$ciWidth

            table$getColumn('dPrimeLower')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))
            table$getColumn('dPrimeUpper')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))
            table$getColumn('cLower')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))
            table$getColumn('cUpper')$setSuperTitle(jmvcore::format('{}% Credible Interval', ciWidth))

            if (is.null(self$options$subj))
                table$addRow(rowKey=1)
        },

        #### Populate tables ----
        .populateSdtGroupTable = function(results) {

            table <- self$results$sdtGroup

            nRows <- nrow(results)
            for (i in 1:nRows)
                table$setRow(rowNo=i, as.list(results[i,]))

        },
        .populateSdtSubjTable = function(results) {

            table <- self$results$sdtSubj

            nRows <- nrow(results)
            for (i in 1:nRows)
                table$setRow(rowNo=i, as.list(results[i,]))

        },

        #### Plot functions ----
        .preparedPrimePlot = function(results) {

            # density <- density(results$samples$d[1,1,,2])
            # df <- data.frame(x = density$x, y = density$y)

            subjs <- split(results$subjs, f = results$subjs$group)
            group <- split(results$group, f = results$group$group)
            typSubj <- split(results$typSubj, f = results$group$group)

            l <- max(sapply(subjs, nrow)):1
            for (i in seq_along(subjs)) {

                typSubj[[i]]$subj <- 'S'
                typSubj[[i]]$x <- l[1] + 1.25
                typSubj[[i]]$subGroup <- 'A'

                subjs[[i]]$x <- l[1:nrow(subjs[[i]])]
                subjs[[i]]$subGroup <- rep('B', nrow(subjs[[i]]))

                subjs[[i]] <- rbind(typSubj[[i]], subjs[[i]][, names(typSubj[[i]])])
            }


            xMin <- min(results$subjs$dPrimeLower, results$typSubj$dPrimeLower, results$group$dPrimeLower)
            xMax <- max(results$subjs$dPrimeUpper, results$typSubj$dPrimeUpper, results$group$dPrimeUpper)

            group <- lapply(group,
                            function(x) {
                                x$xMin = xMin
                                x$xMax = xMax
                                return(x)
                            })

            image <- self$results$dPrimePlot
            image$setState(list(subjs=subjs, group=group))

        },
        .dPrimePlot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            subjs <- image$state$subjs
            group <- image$state$group

            base_breaks_x <- function(breaks, labels) {
                limits <- c(1, max(breaks))
                ggplot2::scale_x_continuous(limits = limits, breaks = breaks, labels = labels)
            }

            base_breaks_y <- function(xMin, xMax, median) {
                values <- pretty(c(xMin, xMax))
                limits <- c(min(values), max(values))
                ggplot2::scale_y_continuous(limits = limits)
            }

            plots <- list()
            for (i in seq_along(subjs)) {

                plots[[i]] <- ggplot2::ggplot(subjs[[i]], ggplot2::aes(x=x, y=dPrime, ymin=dPrimeLower, ymax=dPrimeUpper, color=subGroup)) +
                    ggplot2::geom_pointrange(size = 0.5) + ggplot2::coord_flip() +
                    ggplot2::geom_hline(data = group[[i]], ggplot2::aes(yintercept = 0), linetype = "dotted") +
                    # ggplot2::geom_rect(data = group[[i]], ggplot2::aes(x = NULL, y = NULL, xmin=-Inf, xmax=Inf, ymin=dPrimeLower, ymax=dPrimeUpper, color=NULL), alpha = 0.1) +
                    ggplot2::scale_color_brewer(type="qual", palette = 'Dark2') +
                    ggplot2::scale_fill_brewer(type="qual", palette = 'Dark2') +
                    ggplot2::ggtitle(group[[i]]$group) + ggplot2::theme_bw(base_size = 16) +
                    ggplot2::theme(panel.spacing = grid::unit(2, "lines"), legend.position = "none",
                                   axis.title.y=ggplot2::element_blank(),
                                   axis.title.x=ggplot2::element_blank(),
                                   plot.title = ggplot2::element_text(hjust = 0.5),
                                   panel.border = ggplot2::element_rect(colour = '#333333', fill=NA)) +
                    base_breaks_x(breaks = subjs[[i]]$x, labels = subjs[[i]]$subj) +
                    base_breaks_y(group[[i]]$xMin, group[[i]]$xMax, group[[i]]$dPrime)

            }

            do.call(gridExtra::grid.arrange, c(plots, ncol=length(subjs)))

            # p <- gridExtra::marrangeGrob(plots, ncol=length(subjs), nrow=1)

            # print(p)

            # p <- ggplot2::ggplot(image$state, ggplot2::aes(x=x, y=y)) +
            #     ggplot2::geom_ribbon(ggplot2::aes(ymax=y), ymin=0, fill=theme$fill[2]) +
            #     ggplot2::geom_line(color=theme$color[1]) + ggplot2::xlab('dPrime') + ggplot2::ylab('Density') +
            #     ggtheme



            TRUE

        },

        #### Helper functions ----
        .cleanData = function() {

            df <- self$data
            df <- jmvcore::naOmit(self$data)

            response <- self$options$res
            stimulus <- self$options$stim
            signal <- self$options$sig

            subjects <- self$options$subj
            groups <- self$options$group

            # Calculate number of groups and subjects
            nGroups <- 1
            if ( ! is.null(groups))
                nGroups <- length(unique(df[[groups]]))

            nSubjs <- 1
            if ( ! is.null(subjects)) {
                if (is.null(groups)) {
                    nSubjs <- length(unique(df[[subjects]]))
                } else {
                    nSubjs <- as.numeric(tapply(df[[subjects]],
                                                df[[groups]],
                                                function(x) return(length(unique(x)))))
                }
            }

            # Split data into one data set per subject
            dfs <- df
            if ( ! is.null(subjects)) {
                if ( ! is.null(groups)) {
                    dfs <- split(df, f = list(df[[subjects]], df[[groups]]))
                } else {
                    dfs <- split(df, f = list(df[[subjects]]))
                }
            }

            # Calculate summary stats
            nSignal <- nNoise <- hits <- falseAlarms <- misses <- corrRejects <- matrix(0, max(nSubjs), nGroups)
            iter <- 1
            for (g in 1:nGroups) {
                for (i in 1:nSubjs[g]) {

                    if ( ! is.null(subjects)) {
                        if ( ! is.null(groups)) {
                            respData <- dfs[[paste(iter, g, sep='.')]][[response]]
                            stimData <- dfs[[paste(iter, g, sep='.')]][[stimulus]]
                        } else {
                            respData <- dfs[[iter]][[response]]
                            stimData <- dfs[[iter]][[stimulus]]
                        }
                    } else {
                        respData <- dfs[[response]]
                        stimData <- dfs[[stimulus]]
                    }

                    sigTrials <- which(stimData == signal)
                    noiseTrials <- which(stimData != signal)

                    nSignal[i,g] <- length(sigTrials)
                    nNoise[i,g] <- length(noiseTrials)

                    hits[i,g] <- sum(respData[sigTrials] == stimData[sigTrials])
                    falseAlarms[i,g] <- sum(respData[noiseTrials] != stimData[noiseTrials])
                    misses[i,g] <- nSignal[i,g] - hits[i,g]
                    corrRejects[i,g] <- nNoise[i,g] - falseAlarms[i,g]

                    iter <- iter + 1

                }
            }

            data <- list(nGroups = nGroups, nSubjs = nSubjs, nSignal = nSignal, nNoise = nNoise,
                         HR = hits, FAR = falseAlarms, MR = misses, CR = corrRejects)

            return(data)
        },
        .sumStats = function(data, results) {

            nGroups <- data$nGroups
            nSubjs <- data$nSubjs
            samples <- results$samples

            subjects <- 1
            if ( ! is.null(self$options$subj))
                subjects <- levels(self$data[[self$options$subj]])

            groups <- 1
            if ( ! is.null(self$options$group))
                groups <- levels(self$data[[self$options$group]])

            nRows <- sum(nSubjs)

            group <- subj <- character(nRows)
            dPrime <- dPrimeLower <- dPrimeUpper <- c <- cLower <- cUpper <- h <- fa <- m <- cr <-numeric(nRows)
            dPrimeG <- dPrimeGLower <- dPrimeGUpper <- cG <- cGLower <- cGUpper <- numeric(nGroups)
            dPrimeTyp <- dPrimeTypLower <- dPrimeTypUpper <- cTyp <- cTypLower <- cTypUpper <- numeric(nGroups)

            iter <- 1
            for (g in 1:nGroups) {

                # Group parameters
                dPrimeG[g] <- median(samples$dMu[g,,])
                cG[g] <- median(samples$cMu[g,,])
                cidPrimeG <- HDInterval::hdi(samples$dMu[g,,], self$options$ciWidth / 100)
                ciCG <- HDInterval::hdi(samples$cMu[g,,], self$options$ciWidth / 100)
                dPrimeGLower[g] <- as.numeric(cidPrimeG[1])
                dPrimeGUpper[g] <- as.numeric(cidPrimeG[2])
                cGLower[g] <- as.numeric(ciCG[1])
                cGUpper[g] <- as.numeric(ciCG[2])

                # Typical subject parameters
                dPrimeTyp[g] <- median(samples$dTypSubj[g,,])
                cTyp[g] <- median(samples$cTypSubj[g,,])
                cidPrimeTyp <- HDInterval::hdi(samples$dTypSubj[g,,], self$options$ciWidth / 100)
                ciCTyp <- HDInterval::hdi(samples$cTypSubj[g,,], self$options$ciWidth / 100)
                dPrimeTypLower[g] <- as.numeric(cidPrimeTyp[1])
                dPrimeTypUpper[g] <- as.numeric(cidPrimeTyp[2])
                cTypLower[g] <- as.numeric(ciCTyp[1])
                cTypUpper[g] <- as.numeric(ciCTyp[2])

                # Individual subjects paramters
                for (i in 1:nSubjs[g]) {

                    cidPrime <- HDInterval::hdi(samples$d[i,g,,], self$options$ciWidth / 100)
                    ciC <- HDInterval::hdi(samples$c[i,g,,], self$options$ciWidth / 100)

                    group[iter] <- as.character(groups[g])
                    subj[iter] <- as.character(subjects[iter])

                    dPrime[iter] <- median(samples$d[i,g,,])
                    dPrimeLower[iter] <- as.numeric(cidPrime[1])
                    dPrimeUpper[iter] <- as.numeric(cidPrime[2])
                    c[iter] <- median(samples$c[i,g,,])
                    cLower[iter] <- as.numeric(ciC[1])
                    cUpper[iter] <- as.numeric(ciC[2])

                    h[iter] <- data$HR[i,g]
                    fa[iter] <- data$FAR[i,g]
                    m[iter] <- data$MR[i,g]
                    cr[iter] <- data$CR[i,g]

                    iter <- iter + 1
                }
            }

            dfGroup <- data.frame(group=as.character(groups), dPrime=dPrimeG, dPrimeLower=dPrimeGLower,
                                  dPrimeUpper=dPrimeGUpper, c=cG, cLower=cGLower, cUpper=cGUpper)

            dfTypSubj <- data.frame(group=as.character(groups), dPrime=dPrimeTyp, dPrimeLower=dPrimeTypLower,
                                    dPrimeUpper=dPrimeTypUpper, c=cTyp, cLower=cTypLower, cUpper=cTypUpper)

            dfSubjs <- data.frame(group=group, subj=subj, dPrime=dPrime, dPrimeLower=dPrimeLower,
                                  dPrimeUpper=dPrimeUpper, c=c, cLower=cLower, cUpper=cUpper,
                                  h=h, fa=fa, m=m, cr=cr, stringsAsFactors = FALSE)

            return(list(group=dfGroup, typSubj=dfTypSubj, subjs=dfSubjs))

        })
)
