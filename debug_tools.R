library(R6)

DebugTools <- R6Class(
    public = list(
        use_print_argv = function(argv) {
            message("Start ----\n")
            message("argv <- list()")
            for (arg_name in names(argv)) {
                arg <- argv[[arg_name]]

                if (length(arg) > 1) {
                    arg_string <- paste0("c(\"", paste(arg, collapse="\", \""), "\")")
                }
                else if (is.null(arg)) {
                    next
                    # arg_string <- "NULL"
                }
                else if (is.na(arg)) {
                    arg_string <- "NA"
                }
                else if (is.numeric(arg) || is.logical(arg)) {
                    arg_string <- arg
                }
                else {
                    arg_string <- paste0("\"", arg, "\"")
                }
                message("argv$", arg_name, " <- ", arg_string)
            }
            message("\nEnd ----")
        }
    ),
    private = list()
)

debug_tools <- DebugTools$new()
