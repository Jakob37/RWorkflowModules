library(R6)

DebugTools <- R6Class(
    public = list(
        use_print_argv = function(argv) {
            message("Start ----")
            for (arg_name in names(argv)) {
                arg <- argv[[arg_name]]

                if (length(arg) > 1) {
                    arg_string <- paste0("c(\"", paste(arg, collapse="\", \""), "\")")
                }
                else if (is.null(arg)) {
                    arg_string <- "NULL"
                }
                else if (is.na(arg)) {
                    arg_string <- "NA"
                }
                else {
                    arg_string <- paste0("\"", arg, "\"")
                }
                message("argv$", arg_name, " <- ", arg_string)
            }
            message("End ----")
        }
    ),
    private = list()
)

debug_tools <- DebugTools$new()
