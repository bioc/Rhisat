context("index construction")
test_that("malformed input error",{
    td <- tempdir()
    refs_1 <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),full=TRUE)
    refs_2 <- dir(system.file(package="Rhisat", "extdata", "bt2","reads"),full=TRUE)
    idx <- file.path(td, "lambda_virus")
    # expect_warning(hisat_build(references=NULL, bt2Index=idx))
    expect_error(hisat_build(references=TRUE, bt2Index=idx))
    expect_error(hisat_build(references=refs_1, bt2Index=NULL))
    # expect_warning(hisat_build(references=refs_2, bt2Index=idx))
    # expect_warning(hisat_build(references=refs_1, bt2Index=idx,"--noparam"))
}
)

test_that("formatted input work",{
    td <- tempdir()
    idx <- file.path(td, "hisat-build-test")
    refs <- dir(system.file(package="Rhisat", "extdata", "bt2","refs"),full=TRUE)
    expect_equal(
        hisat_build(references=refs, bt2Index=idx, exe=FALSE),
        paste(
            file.path(system.file(package="Rhisat"), "hisat-build"),
            refs, idx
        )
    )
}
)