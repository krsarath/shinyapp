library(testthat)

app_dir <- ".."
repo_dir <- "../.."

test_that("app.R file exists in repository", {
  expect_true(file.exists(file.path(app_dir, "app.R")))
})

test_that("Required packages load without error", {
  expect_no_error(library(shiny))
  expect_no_error(library(DESeq2))
  expect_no_error(library(edgeR))
  expect_no_error(library(NOISeq))
  expect_no_error(library(ggplot2))
  expect_no_error(library(ggVennDiagram))
})

test_that("LICENSE file exists", {
  expect_true(file.exists(file.path(repo_dir, "LICENSE")))
})

test_that("CONTRIBUTING.md exists", {
  expect_true(file.exists(file.path(repo_dir, "CONTRIBUTING.md")))
})

test_that("www folder exists for app assets", {
  expect_true(dir.exists(file.path(app_dir, "www")))
})
