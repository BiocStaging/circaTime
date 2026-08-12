# GitHub-only engines (tauFisher, TimeSignatR, zeitzeiger) are `Suggests`, not
# `Imports`: Bioconductor/CRAN flag GitHub-only strong dependencies (Remotes on
# Imports/Depends) as a submission blocker. Declaring them as Suggests keeps the
# other engines (Molecular Timetable) usable with a base install; each wrapper
# checks for its own package at call time and fails with install instructions
# rather than a cryptic "could not find function" error.

.requireEngine <- function(pkg, install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      pkg, " is required for this engine but is not installed.\n",
      "Install it with: ", install,
      call. = FALSE
    )
  }
}
