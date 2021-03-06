
# This file is automatically generated, you probably don't want to edit this

resWagnerOptions <- if (requireNamespace('jmvcore')) R6::R6Class(
    "resWagnerOptions",
    inherit = jmvcore::Options,
    public = list(
        initialize = function(
            choice = NULL,
            key = NULL, ...) {

            super$initialize(
                package='resWagner',
                name='resWagner',
                requiresData=TRUE,
                ...)

            private$..choice <- jmvcore::OptionVariable$new(
                "choice",
                choice)
            private$..key <- jmvcore::OptionVariable$new(
                "key",
                key)

            self$.addOption(private$..choice)
            self$.addOption(private$..key)
        }),
    active = list(
        choice = function() private$..choice$value,
        key = function() private$..key$value),
    private = list(
        ..choice = NA,
        ..key = NA)
)

resWagnerResults <- if (requireNamespace('jmvcore')) R6::R6Class(
    inherit = jmvcore::Group,
    active = list(
        plot = function() private$.items[["plot"]]),
    private = list(),
    public=list(
        initialize=function(options) {
            super$initialize(
                options=options,
                name="",
                title="Rescorla-Wagner Model")
            self$add(jmvcore::Image$new(
                options=options,
                name="plot",
                title="",
                renderFun=".plot",
                clearWith=list(
                    "choice",
                    "key")))}))

resWagnerBase <- if (requireNamespace('jmvcore')) R6::R6Class(
    "resWagnerBase",
    inherit = jmvcore::Analysis,
    public = list(
        initialize = function(options, data=NULL, datasetId="", analysisId="", revision=0) {
            super$initialize(
                package = 'resWagner',
                name = 'resWagner',
                version = c(1,0,0),
                options = options,
                results = resWagnerResults$new(options=options),
                data = data,
                datasetId = datasetId,
                analysisId = analysisId,
                revision = revision,
                pause = NULL,
                completeWhenFilled = FALSE)
        }))

#' Rescorla-Wagner Model
#'
#' 
#' @param data .
#' @param choice .
#' @param key .
#' @return A results object containing:
#' \tabular{llllll}{
#'   \code{results$plot} \tab \tab \tab \tab \tab a density plot \cr
#' }
#'
#' @export
resWagner <- function(
    data,
    choice,
    key) {

    if ( ! requireNamespace('jmvcore'))
        stop('resWagner requires jmvcore to be installed (restart may be required)')

    options <- resWagnerOptions$new(
        choice = choice,
        key = key)

    results <- resWagnerResults$new(
        options = options)

    analysis <- resWagnerClass$new(
        options = options,
        data = data)

    analysis$run()

    analysis$results
}
