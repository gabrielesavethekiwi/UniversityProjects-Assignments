PROGRAM NUMBER_PRECISION
implicit none


INTEGER*2 :: small_2, big_2 
INTEGER*4 :: small_4, big_4

REAL*4 :: PI_4	!IEEE standard: single precision = 32 bits = 4 bytes
REAL*8 :: PI_8  !double prec = 64 bits


small_2 = 1
big_2 = 2000000 !overflow, add -fno-range-check to compiler command

print *, '2 bytes: ', small_2+big_2

!repeat procedure with 4 bytes 
small_4 = 1
big_4 = 2000000

print *, '4 bytes: ', small_4+big_4   !no overflow

PI_4 = 4.0_4 * atan(1.0_4)
PI_8 = 4.0_8 * atan(1.0_8)


print *, 'single precision: ', 1e32*PI_4 + sqrt(2.0_4)*1e21   
print *, 'double precision: ', 1e32*PI_8 + sqrt(2.0_8)*1e21   

!the first sum doesn't change pi*10^32 due to lack of precision, while from the second we recover the first 5 digits of sqrt2 correctly

print *, 'second addendum invisible: ', (1e32*PI_4 + sqrt(2.0_4)*1e21) - 1e32*PI_4, sqrt(2.0_4)*1e21 
print *, 'second addendum detected (5 digits):', (1e32*PI_8 + sqrt(2.0_8)*1e21) - 1e32*PI_8, sqrt(2.0_8)*1e21

END PROGRAM number_precision