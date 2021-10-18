vec1 <- c(1, 0)
vec2 <- c(0,1)
circle <- .construct_circle(vec1, vec2)
percentage_val <- 0
enforce_boundary = T

rad1 <- .find_radian(circle, vec1)
rad2 <- .find_radian(circle, vec2)
stopifnot(rad1 < 0 & rad2 > 0) # must be true based on how we constructed vec1 and vec2
rad1 <- rad1 + 2*pi # to ensure rad1 is larger than rad2

rad1*180/pi
360-45
rad2*180/pi
90+45

stopifnot(abs(sum(vec1 - c(1,0))) <= 1e-6,
          circle$radius >= circle$center - 1e-6) # must be true based on how we constructed vec1
position1 <- circle$center[1] - sqrt(circle$radius^2 - circle$center[2]^2) 
rad1_truncated <- .find_radian(circle, c(position1, 0))
if(rad1_truncated < 0) rad1_truncated <- rad1_truncated + 2*pi
rad1_truncated*180/pi
180+45

inner_rad <- atan(circle$center[1]/circle$center[2])
inner_rad*180/pi
if(inner_rad < 0) inner_rad <- inner_rad + 2*pi
radmid_truncated <- 3*pi/2 - inner_rad
radmid_truncated*180/pi
stopifnot(rad1_truncated > radmid_truncated)
rad2_truncated <- radmid_truncated - (rad1_truncated - radmid_truncated)
rad2_truncated*180/pi
