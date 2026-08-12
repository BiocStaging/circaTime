# S4 generics/methods (AGENTS.md Section 5/10). `show` is already a generic from
# the methods package; only the method needs defining here.

#' @rdname PhaseReferenceModel-class
#' @param object A `PhaseReferenceModel`.
#' @importFrom methods show
#' @export
methods::setMethod("show", "PhaseReferenceModel", function(object) {
  cat(sprintf(
    "PhaseReferenceModel: species=%s tissue=%s method=%s platform=%s\n",
    object@species, object@tissue, object@method, object@platform
  ))
  cat(sprintf(
    "  %d amplitude gene(s), %d residual-pool draw(s), ampThreshold=%.3f\n",
    length(object@ampGenes), length(object@residualPool), object@ampThreshold
  ))
  cat(sprintf(
    "  held-out CV circular MAE: %s\n",
    if (is.na(object@cvMAE)) "not measured" else sprintf("%.2fh", object@cvMAE)
  ))
  cat(sprintf("  trained on: %s\n", object@trainedOn))
  invisible(object)
})
