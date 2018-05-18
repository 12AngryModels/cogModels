
resWagnerClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "resWagnerClass",
    inherit = resWagnerBase,
    private = list(
        .run = function() {

            path <- system.file("models/resWagner.txt", package = "resWagner")
            
            choice <- self$options$choice
            key <- self$options$key
            
            if (is.null(choice) || is.null(key))
                return()
            
            choiceData <- as.numeric(self$data[[choice]])
            keyData <- as.numeric(self$data[[key]])
            rewardData <- as.numeric(choiceData == keyData)
            
            nTrials <- length(choiceData)
            
            data <- list('choice'=choiceData, 'reward'=rewardData, 'N'=nTrials)
            
            model <- rjags::jags.model(file = path,
                                       data = data,
                                       n.chains = 4,
                                       n.adapt = 100)
            
            update(model, 1000)
            
            samples <- rjags::jags.samples(model, c('a'), 1000)
            
            image <- self$results$plot
            image$setState(samples[['a']])
            
        },
        .plot = function(image, ggtheme, theme, ...) {
            
            if (is.null(image$state))
                return(FALSE)
            
            plot(density(image$state))
            
            TRUE
    
        }
))
