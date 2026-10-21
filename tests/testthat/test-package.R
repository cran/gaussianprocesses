test_that("package scaffold is loadable", {
  expect_s3_class(
    packageVersion("gaussianprocesses"),
    "package_version"
  )
})


test_that("exported defaults do not reference internal helpers", {
  namespace <- asNamespace("gaussianprocesses")
  exports <- getNamespaceExports(namespace)
  internal <- setdiff(
    ls(namespace, all.names = TRUE),
    exports
  )

  for (name in exports) {
    used <- as.character(
      unlist(
        lapply(
          formals(get(name, envir = namespace)),
          all.names
        ),
        use.names = FALSE
      )
    )

    expect_identical(
      intersect(used, internal),
      character(0),
      info = name
    )
  }
})


test_that("public numerical tolerances equal the internal default", {
  tolerance_arguments <- list(
    fit_gp = "symmetry_tolerance",
    optimize_gp = "symmetry_tolerance",
    fit_sparse_gp = "symmetry_tolerance",
    simulate_gp_data = "symmetry_tolerance",
    sample_gp_prior = "symmetry_tolerance",
    fit_heteroscedastic_gp = "symmetry_tolerance",
    predict_gp = "variance_tolerance",
    predict_sparse_gp = "variance_tolerance",
    predict_heteroscedastic_gp = "variance_tolerance",
    forecast_gp = "variance_tolerance"
  )

  for (name in names(tolerance_arguments)) {
    argument <- tolerance_arguments[[name]]
    default <- formals(
      getExportedValue("gaussianprocesses", name)
    )[[argument]]

    expect_identical(
      eval(default, baseenv()),
      gaussianprocesses:::.default_tolerance(),
      info = paste0(name, "(", argument, ")")
    )
  }
})
