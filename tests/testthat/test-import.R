test_that("import_examples_pH_data works", {
    
    
    read_pH(
        system.file( "extdata", "Igrashi_copy.csv", package = "GatorpH")
    ) |> expect_no_error()

    system.file( "extdata", "Birkhed1_copy.csv", package = "GatorpH") |> 
        read_pH() |> 
        expect_no_error()

})


test_that("import_template_data works", {
    
    dumdat1 <- system.file("extdata", "template.csv", package = "GatorpH") |> 
        read.csv() |> 
        rbind(data.frame(
            id = c(1, 1, 1, 1, 1, 1, 1),
            time = c(1, 5, 10, 15, 20, 30, 50),
            pH = c(7, 6.5, 6.2, 6.0, 5.8, 5.5, 5.2),
            group = "A",
            baseline = 7
        )) 

    dumdat2 <- dumdat1
    dumdat2$baseline <- NA

    # write to temp file
    tmpfile1 <- tempfile(fileext = ".csv")
    tmpfile2 <- tempfile(fileext = ".csv")

    write.csv(dumdat1, tmpfile1, row.names = FALSE)
    write.csv(dumdat2, tmpfile2, row.names = FALSE)
    # test import
    read_pH(tmpfile1) |> expect_no_error()
    read_pH(tmpfile2) |> expect_no_error()

})
