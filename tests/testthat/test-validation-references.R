# Independent references for the inference methods ----------------------------
#
# inst/validation/references.R defines the references; the validation
# article tabulates them. Each reference is one test here, skipped when its
# reference package is not installed.

source(
  system.file("validation", "references.R", package = "gaussianprocesses"),
  local = TRUE
)


test_that("every exported inference function belongs to a method with a reference", {
  exports <- getNamespaceExports("gaussianprocesses")
  covered <- unlist(validation_methods, use.names = FALSE)
  references <- validation_references()
  methods <- vapply(references, function(reference) reference$method, "")

  expect_identical(
    sort(setdiff(exports, c(covered, names(validation_exclusions)))),
    character(0)
  )
  expect_identical(
    setdiff(c(covered, names(validation_exclusions)), exports),
    character(0)
  )
  expect_false(anyDuplicated(covered) > 0)
  expect_setequal(unique(methods), names(validation_methods))
})


test_that("references are well formed", {
  references <- validation_references()
  ids <- vapply(references, function(reference) reference$id, "")
  suggested <- trimws(sub(
    "\\(.*", "",
    strsplit(utils::packageDescription("gaussianprocesses")$Suggests,
             ",")[[1]]
  ))

  expect_false(anyDuplicated(ids) > 0)

  for (reference in references) {
    expect_true(
      reference$type %in% c("implementation", "brute force", "closed form"),
      label = reference$id
    )
    expect_true(is.numeric(reference$tolerance) && reference$tolerance >= 0,
                label = reference$id)
    expect_true(is.function(reference$run), label = reference$id)
    expect_true(nzchar(reference$quantity) && nzchar(reference$reference),
                label = reference$id)

    if (!is.na(reference$package)) {
      expect_true(reference$package %in% suggested, label = reference$id)
      expect_false(is.na(reference$recorded_version), label = reference$id)
    }
  }
})


for (reference in validation_references()) {
  local({
    reference <- reference

    test_that(sprintf("reference %s agrees within its tolerance", reference$id), {
      # The references take a few seconds in all; they run everywhere
      # except on CRAN, to keep its check within about a minute.
      skip_on_cran()

      if (!is.na(reference$package)) {
        skip_if_not_installed(reference$package)
      }

      expect_lte(reference$run(), reference$tolerance, label = reference$id)
    })
  })
}


test_that("references whose package is missing are reported, not run", {
  missing <- validation_reference(
    "missing-package", "Exact regression", "a quantity", "a package",
    "implementation", tolerance = 1, package = "aPackageThatIsNotInstalled",
    recorded_version = "1.0", run = function() stop("must not run")
  )
  result <- run_validation_references(list(missing))

  expect_identical(result$status, "not run")
  expect_true(is.na(result$discrepancy))
  expect_true(is.na(result$installed_version))
})
