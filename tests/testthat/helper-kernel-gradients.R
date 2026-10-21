# Helpers for checking analytical kernel derivatives against central finite
# differences in log-parameter space.

# Step in log-parameter space and the tolerance on the scaled discrepancy
# max|analytical - numerical| / max(1, max|analytical|). With this step the
# central-difference truncation error is about 1e-10 for the cases below.
gradient_fd_step <- 1e-6
gradient_fd_tolerance <- 1e-7


finite_difference_kernel_gradient <- function(
    kernel,
    x,
    y = NULL,
    step = gradient_fd_step) {
  parameters <- kernel_parameters(kernel, flatten = TRUE)

  derivatives <- lapply(
    seq_along(parameters),
    function(index) {
      covariance_at <- function(direction) {
        shifted <- parameters[index] * exp(direction * step)
        evaluate_kernel(update_kernel_parameters(kernel, shifted), x, y)
      }

      (covariance_at(1) - covariance_at(-1)) / (2 * step)
    }
  )
  names(derivatives) <- names(parameters)
  derivatives
}


scaled_gradient_discrepancy <- function(analytical, numerical) {
  max(abs(analytical - numerical)) / max(1, max(abs(analytical)))
}


# Representative kernels, parameter values, and inputs: 1D, multidimensional
# with shared and ARD length scales, cross-covariances, coincident inputs,
# and nested composites.
kernel_gradient_cases <- function() {
  x_1d <- c(-1.3, -0.4, 0, 0.35, 1.2, 2.1)
  y_1d <- c(-0.8, 0.5, 1.7, 3)
  x_coincident <- c(0, 0, 0.5, 1, 1)
  x_2d <- cbind(seq(-1, 1, length.out = 7), cos(1:7))
  y_2d <- cbind(c(0.2, -0.7, 1.4, 0.1, 2), sin(1:5))
  x_3d <- cbind(
    seq(-1, 1, length.out = 6),
    sin(1:6),
    c(0.3, -1.1, 0.8, 0, 2.2, -0.4)
  )
  y_3d <- x_3d[c(2, 5, 6), ] + 0.25

  case <- function(name, kernel, x, y = NULL) {
    list(name = name, kernel = kernel, x = x, y = y)
  }

  stationary <- list(
    rbf = rbf_kernel,
    matern12 = matern12_kernel,
    matern32 = matern32_kernel,
    matern52 = matern52_kernel
  )

  cases <- list()

  for (type in names(stationary)) {
    constructor <- stationary[[type]]
    cases <- c(
      cases,
      list(
        case(paste(type, "1D"), constructor(1, 0.7), x_1d),
        case(paste(type, "1D cross"), constructor(2.5, 3), x_1d, y_1d),
        case(paste(type, "coincident"), constructor(0.6, 0.4), x_coincident),
        case(
          paste(type, "2D ARD"),
          constructor(0.8, c(0.4, 2)),
          x_2d
        ),
        case(paste(type, "3D shared"), constructor(1.3, 1.1), x_3d),
        case(
          paste(type, "3D ARD cross"),
          constructor(1.7, c(0.5, 1.5, 3)),
          x_3d,
          y_3d
        )
      )
    )
  }

  c(
    cases,
    list(
      case(
        "rational quadratic small alpha",
        rational_quadratic_kernel(1.2, 0.8, alpha = 0.3),
        x_1d
      ),
      case(
        "rational quadratic moderate alpha cross",
        rational_quadratic_kernel(0.5, 1.6, alpha = 2),
        x_1d,
        y_1d
      ),
      case(
        "rational quadratic large alpha",
        rational_quadratic_kernel(1, 0.9, alpha = 50),
        x_1d
      ),
      case(
        "rational quadratic 2D ARD",
        rational_quadratic_kernel(0.9, c(0.3, 1.8), alpha = 1.5),
        x_2d
      ),
      case(
        "rational quadratic 3D ARD cross",
        rational_quadratic_kernel(1.4, c(2, 0.6, 1), alpha = 0.8),
        x_3d,
        y_3d
      ),
      case("periodic 1D", periodic_kernel(1, 0.8, period = 1), x_1d),
      case(
        "periodic 1D cross",
        periodic_kernel(2, 1.7, period = 2.5),
        x_1d,
        y_1d
      ),
      case(
        "periodic 2D ARD",
        periodic_kernel(0.8, c(0.6, 2), period = 1.3),
        x_2d
      ),
      case(
        "periodic 3D shared cross",
        periodic_kernel(1.1, 0.9, period = 0.7),
        x_3d,
        y_3d
      ),
      case("linear 1D", linear_kernel(0.4), x_1d),
      case("linear 3D cross", linear_kernel(2.3), x_3d, y_3d),
      case("white noise coincident", white_noise_kernel(0.2), x_coincident),
      case(
        "white noise cross",
        white_noise_kernel(1.5),
        x_coincident,
        c(0, 1, 2)
      ),
      case(
        "sum ARD plus periodic",
        sum_kernel(
          rbf_kernel(0.8, c(0.4, 2)),
          periodic_kernel(0.5, 1.2, period = 1.1)
        ),
        x_2d
      ),
      case(
        "product rbf linear",
        product_kernel(rbf_kernel(1.2, 0.9), linear_kernel(0.4)),
        x_1d
      ),
      case(
        "product of three cross",
        product_kernel(
          matern32_kernel(1, 1.3),
          periodic_kernel(0.7, 0.8, period = 2),
          linear_kernel(0.5)
        ),
        x_1d,
        y_1d
      ),
      case(
        "product with white noise zeros",
        product_kernel(rbf_kernel(1.1, 0.6), white_noise_kernel(0.3)),
        x_coincident
      ),
      case(
        "scaled ARD matern52",
        scale_kernel(matern52_kernel(1.4, c(0.7, 1.9)), scale = 2.5),
        x_2d,
        y_2d
      ),
      case(
        "nested scale sum product",
        scale_kernel(
          sum_kernel(
            product_kernel(
              rbf_kernel(0.9, c(0.5, 1.5)),
              periodic_kernel(1.2, 0.7, period = 1.6)
            ),
            white_noise_kernel(0.05)
          ),
          scale = 0.7
        ),
        x_2d
      ),
      case(
        "sum of scaled kernels 3D ARD",
        sum_kernel(
          scale_kernel(rbf_kernel(1, c(0.8, 1.2, 2)), scale = 1.3),
          scale_kernel(
            rational_quadratic_kernel(0.6, 0.9, alpha = 0.7),
            scale = 0.4
          )
        ),
        x_3d,
        y_3d
      ),
      case(
        "time-series kernel",
        time_series_kernel(
          period = 2,
          trend_variance = 0.05,
          periodic_variance = 1,
          periodic_length_scale = 0.8,
          local_variance = 0.3,
          local_length_scale = 1.5
        ),
        x_1d
      )
    )
  )
}


# Central differences of the covariance in each parameter's optimizer
# coordinate, which is what kernel_gradient() differentiates with respect to.
coordinate_fd_kernel_gradient <- function(kernel, x, y = NULL, step = 1e-6) {
  parameters <- kernel_parameters(kernel, flatten = TRUE)
  constraints <- gaussianprocesses:::.flatten_kernel_constraints(kernel)

  derivatives <- lapply(
    names(parameters),
    function(path) {
      constraint <- constraints[[path]]
      coordinate <- constraint$to_coordinate(parameters[[path]])

      covariance_at <- function(direction) {
        shifted <- parameters[path]
        shifted[] <- constraint$from_coordinate(coordinate + direction * step)
        evaluate_kernel(update_kernel_parameters(kernel, shifted), x, y)
      }

      (covariance_at(1) - covariance_at(-1)) / (2 * step)
    }
  )
  names(derivatives) <- names(parameters)
  derivatives
}
