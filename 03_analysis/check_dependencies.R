required_packages <- c(
  survival = "3.8-6",
  ggplot2 = "4.0.3",
  patchwork = "1.3.2"
)

missing <- names(required_packages)[
  !vapply(names(required_packages), requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop(
    "Missing required package(s): ",
    paste(missing, collapse = ", "),
    ". Install them before running the pipeline."
  )
}

for (pkg in names(required_packages)) {
  installed <- as.character(packageVersion(pkg))
  message(pkg, ": installed ", installed, "; verified analysis version ", required_packages[[pkg]])
}

message("R: installed ", getRversion(), "; verified analysis version 4.6.1")
