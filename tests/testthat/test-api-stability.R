# API stability tiers and deprecations (R/stability.R) -------------------------

exports <- getNamespaceExports("gaussianprocesses")

# Help pages as Rd source text, from the source tree when the tests run there
# and from the installed package otherwise, as under R CMD check.
help_pages <- function() {
  man <- test_path("..", "..", "man")
  db <- if (dir.exists(man)) {
    files <- list.files(man, "\\.Rd$", full.names = TRUE)
    stats::setNames(lapply(files, tools::parse_Rd), basename(files))
  } else {
    tools::Rd_db("gaussianprocesses")
  }
  names(db) <- sub("\\.Rd$", "", names(db))
  vapply(db, function(rd) {
    paste(as.character(rd, deparse = TRUE), collapse = "")
  }, character(1))
}

page_aliases <- function(text) {
  regmatches(text, gregexpr("(?<=\\\\alias\\{)[^}]+", text, perl = TRUE))[[1]]
}

page_tier <- function(text) {
  tier <- regmatches(
    text,
    regexpr("(?<=\\\\section\\{Stability\\}\\{)\\s*[A-Z][a-z]+", text,
            perl = TRUE)
  )
  if (length(tier) == 0L) NA_character_ else tolower(trimws(tier))
}


test_that("the tier table names exports only", {
  experimental <- gaussianprocesses:::.experimental_exports
  deprecated <- names(gaussianprocesses:::.deprecated_exports)

  expect_false(anyDuplicated(experimental) > 0)
  expect_identical(setdiff(c(experimental, deprecated), exports), character(0))
  expect_identical(intersect(experimental, deprecated), character(0))
  expect_true(all(gaussianprocesses:::.deprecated_exports %in% exports))
})


test_that("every export's help page states its tier", {
  pages <- help_pages()
  aliases <- lapply(pages, page_aliases)

  for (name in exports) {
    topic <- names(aliases)[vapply(aliases, function(a) name %in% a, TRUE)]

    expect_length(topic, 1L)
    expect_identical(page_tier(pages[[topic]]),
                     gaussianprocesses:::.api_tier(name),
                     label = name)
  }
})


test_that("the reference index is grouped by tier", {
  config <- test_path("..", "..", "_pkgdown.yml")
  skip_if_not(file.exists(config), "_pkgdown.yml is not in the package")
  skip_if_not_installed("yaml")

  pages <- help_pages()
  aliases <- lapply(pages, page_aliases)
  tiers <- c("Stable interfaces" = "stable",
             "Experimental interfaces" = "experimental",
             "Deprecated" = "deprecated")
  current <- NA_character_
  listed <- character(0)

  for (entry in yaml::read_yaml(config)$reference) {
    if (!is.null(entry$title)) {
      expect_true(entry$title %in% names(tiers), label = entry$title)
      current <- tiers[[entry$title]]
    }

    for (topic in entry$contents) {
      documented <- intersect(aliases[[topic]], exports)
      expect_true(topic %in% names(pages), label = topic)
      expect_true(
        all(gaussianprocesses:::.api_tier(documented) == current),
        label = paste(topic, "under", current)
      )
      listed <- c(listed, topic)
    }
  }

  with_exports <- names(aliases)[vapply(aliases, function(a) {
    any(a %in% exports)
  }, TRUE)]
  expect_false(anyDuplicated(listed) > 0)
  expect_identical(setdiff(with_exports, listed), character(0))
})


test_that("the deprecation tables match the code and the help page", {
  arguments <- gaussianprocesses:::.deprecated_arguments
  fields <- gaussianprocesses:::.deprecated_fields
  deprecated_exports <- gaussianprocesses:::.deprecated_exports
  namespace <- asNamespace("gaussianprocesses")

  for (i in seq_len(nrow(arguments))) {
    expect_true(
      arguments$argument[[i]] %in%
        names(formals(get(arguments$fun[[i]], envir = namespace))),
      label = paste(arguments$fun[[i]], arguments$argument[[i]])
    )
  }

  expect_named(arguments, c("fun", "argument", "replacement"))
  expect_named(fields, c("fun", "field", "replacement"))
  expect_true(all(deprecated_exports %in% exports))

  # A deprecated name is documented on ?gaussianprocesses-deprecated, which
  # exists only while there are any (none in 1.0.0).
  names <- c(names(deprecated_exports), arguments$argument, fields$field)

  if (length(names) > 0L) {
    deprecated_page <- help_pages()[["gaussianprocesses-deprecated"]]

    for (name in c(names, deprecated_exports, arguments$replacement,
                   fields$replacement)) {
      expect_match(deprecated_page, name, fixed = TRUE, label = name)
    }
  }
})


test_that("names removed in 1.0.0 do not return", {
  removed <- gaussianprocesses:::.removed_in_1.0.0
  namespace <- asNamespace("gaussianprocesses")

  expect_false(any(names(removed$exports) %in% exports))
  expect_false(any(names(removed$exports) %in% ls(namespace)))

  functions <- unique(c(
    exports,
    grep("^(sample_gp_posterior|predict)\\.", ls(namespace), value = TRUE)
  ))
  for (name in functions) {
    value <- get(name, envir = namespace)

    if (is.function(value)) {
      expect_length(
        intersect(names(formals(value)), names(removed$arguments)),
        0L
      )
    }
  }

  expect_error(
    fit_gp(seq(0, 1, length.out = 5), 1:5, rbf_kernel(), tolerance = 1e-6),
    "unused argument"
  )

  simulation <- gp_simulation_scenario("heteroscedastic_1d", n = 20)
  model <- fit_heteroscedastic_gp(
    simulation$x, simulation$observed,
    kernel = matern52_kernel(length_scale = 0.2),
    noise_kernel = rbf_kernel(length_scale = 0.4),
    max_iterations = 30
  )
  prediction <- predict_heteroscedastic_gp(model, c(0.2, 0.8))
  expect_false(any(names(removed$fields) %in% names(prediction)))
  expect_false(is.null(prediction[[removed$fields[["noise_variance"]]]]))
})
