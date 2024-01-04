test_that("rotation of coordinate system works", {
  pt <- matrix(c(1, 0), nrow = 2)
  
  ############################
  # right-handed, clock-wise #
  ############################
  
  # 30 degrees
  expect_equal(
    get_rot_mtx(pi / 6, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(3) / 2, 0.5), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 6, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(3) / 2, 0.5), nrow = 2)
  )
  
  # 45 degrees
  expect_equal(
    get_rot_mtx(pi / 4, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(2) / 2, sqrt(2) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 4, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(2) / 2, sqrt(2) / 2), nrow = 2)
  )
  
  # 60 degrees
  expect_equal(
    get_rot_mtx(pi / 3, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(0.5, sqrt(3) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 3, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(0.5, sqrt(3) / 2), nrow = 2)
  )
  
  # 90 degrees
  expect_equal(
    get_rot_mtx(pi / 2, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(0, 1), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 2, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(0, 1), nrow = 2)
  )
  
  # 180 degrees
  expect_equal(
    get_rot_mtx(pi, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
  
  ###################################
  # right-handed, counterclock-wise #
  ###################################
  
  # 30 degrees
  expect_equal(
    get_rot_mtx(-pi / 6, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(3) / 2, -0.5), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 6, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(3) / 2, -0.5), nrow = 2)
  )
  
  # 45 degrees
  expect_equal(
    get_rot_mtx(-pi / 4, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(2) / 2, -sqrt(2) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 4, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(2) / 2, -sqrt(2) / 2), nrow = 2)
  )

  # 60 degrees
  expect_equal(
    get_rot_mtx(-pi / 3, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(0.5, -sqrt(3) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 3, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(0.5, -sqrt(3) / 2), nrow = 2)
  )

  # 90 degrees
  expect_equal(
    get_rot_mtx(-pi / 2, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(0, -1), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 2, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(0, -1), nrow = 2)
  )

  # 180 degrees
  expect_equal(
    get_rot_mtx(-pi, is.rh = TRUE, is.clock = TRUE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi, is.rh = TRUE, is.clock = FALSE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
  
  ###########################
  # left-handed, clock-wise #
  ###########################
  
  # 30 degrees
  expect_equal(
    get_rot_mtx(pi / 6, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(3) / 2, -0.5), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 6, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(3) / 2, -0.5), nrow = 2)
  )
  
  # 45 degrees
  expect_equal(
    get_rot_mtx(pi / 4, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(2) / 2, -sqrt(2) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 4, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(2) / 2, -sqrt(2) / 2), nrow = 2)
  )

  # 60 degrees
  expect_equal(
    get_rot_mtx(pi / 3, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(0.5, -sqrt(3) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 3, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(0.5, -sqrt(3) / 2), nrow = 2)
  )

  # 90 degrees
  expect_equal(
    get_rot_mtx(pi / 2, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(0, -1), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi / 2, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(0, -1), nrow = 2)
  )

  # 180 degrees
  expect_equal(
    get_rot_mtx(pi, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(-pi, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
  
  ##################################
  # left-handed, counterclock-wise #
  ##################################
  
  # 30 degrees
  expect_equal(
    get_rot_mtx(-pi / 6, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(3) / 2, 0.5), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 6, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(3) / 2, 0.5), nrow = 2)
  )
  
  # 45 degrees
  expect_equal(
    get_rot_mtx(-pi / 4, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(sqrt(2) / 2, sqrt(2) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 4, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(sqrt(2) / 2, sqrt(2) / 2), nrow = 2)
  )

  # 60 degrees
  expect_equal(
    get_rot_mtx(-pi / 3, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(0.5, sqrt(3) / 2), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 3, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(0.5, sqrt(3) / 2), nrow = 2)
  )

  # 90 degrees
  expect_equal(
    get_rot_mtx(-pi / 2, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(0, 1), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi / 2, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(0, 1), nrow = 2)
  )

  # 180 degrees
  expect_equal(
    get_rot_mtx(-pi, is.rh = FALSE, is.clock = TRUE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
  expect_equal(
    get_rot_mtx(pi, is.rh = FALSE, is.clock = FALSE) %*% pt,
    matrix(c(-1, 0), nrow = 2)
  )
})
