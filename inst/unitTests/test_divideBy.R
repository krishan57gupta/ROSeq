test_divideBy <- function() {
    checkEquals(divideBy(8, 2), 2)
    checkTrue(is.na(divideBy(10, 0)))
    checkEqualsNumeric(divideBy(5, 1.2345), 3.24, tolerance=1.0e-4)
}