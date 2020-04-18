test_divideBy <- function() {
    checkEquals(divideBy(8, 4), 2)
    checkTrue(is.na(divideBy(10, 0)))
}
