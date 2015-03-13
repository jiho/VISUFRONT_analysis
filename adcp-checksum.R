x <- "GPHDT,-999.990,T" # check = 11
x <- "GPHDT,999.000,T"  # check = 3C
x <- "GPHDT,201.159,T"  # check = 3B
x <- "GPHDT,199.867,T"   # check = 3D

# checksum for end of lines in NMEA records
nmea_check <- function(x) {
  r <- charToRaw(x)
  check <- r[1]
  for (i in 2:length(r)) {
    check <- xor(check, r[i])
  }
  return(check)
}

# OK, works