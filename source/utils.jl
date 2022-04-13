#############################################################################
#### FRAME CHANGING UTILITIES
#############################################################################
"""
    sph2cart(r,theta,phi)

Spherical to cartesian coordinate transform

# Arguments
- `r::Float` : spherical radius
- `theta::Float` : longitude [0,2*pi]
- `phi::Float` : colatitude [0,pi]

# Output
- `x,y,z`: cartesian coordinates (same units as r)

# Dependency
None

"""
function sph2cart(r,theta,phi)
    x = r* sin(phi)* cos(theta)
    y = r* sin(phi)* sin(theta)
    z = r* cos(phi)
    return x, y, z
end


"""
    cart2sph(x, y, z)

Cartesian to spherical transform

# Arguments
- `x::Number` : x coordinate
- `y::Number` : y coordinate
- `z::Number` : z coordinate

# Output
- `r ,theta, phi` : spherical coordinates, spherical radius, longitude [0,2*pi], colatitude [0,pi]

# Dependency
None 

"""
function cart2sph(x, y, z)
    r = sqrt(x*x + y*y + z*z)
    theta = atan(y, x)
    if theta < 0
        theta = theta + 2*pi
    end
    phi = acos(z/r)
    return r, theta, phi
end


"""
    cart2cyl(x,y,z)

Cartesian to cylindrical coordinate transform

# Arguments
- `x::Number` : x coordinate
- `y::Number` : y coordinate
- `z::Number` : z coordinate

# Output
- `R, theta, z` : cylindrical coordinates, cylindrical radius R, longitude theta, altitude z

# Dependency
None

"""
function cart2cyl(x, y, z)
    R = sqrt(x*x + y*y)
    theta = atan(y, x)
    if theta<0
        theta = theta + 2*pi
    end
    return R, theta, z
end

"""
    cyl2cart(R, theta, z)

Cylindrical to cartesian transform

# Arguments
- `R::Number` : cylindrical radius, distance in the equatorial plane
- `theta::Number` : longitude, angle in the equatorial plane [0., 2*pi]
- `z::Number` : altitude, z coordinate

# Output
- `x, y, z` : cartesian coordinates

# Dependency
None
"""
function cyl2cart(R, theta, z)
    x = R * cos(theta)
    y = R * sin(theta)
    return x, y, z
end

"""
    cyl2sph(R, theta, z)

Cylindrical to cartesian transform

# Arguments
- `R::Number` : cylindrical radius, distance in the equatorial plane
- `theta::Number` : longitude, angle in the equatorial plane [0., 2*pi]
- `z::Number` : altitude, z coordinate

# Output
- `r, theta, phi` : spherical coordinates, r spherical radius, theta longitude [0, 2*pi], phi colatitude [0,pi]

# Dependency
None
"""
function cyl2sph(R, theta, z)
    r = sqrt(R*R+z*z)
    phi = pi/2. - atan(z,R)
    return r, theta, phi
end

"""
    sph2cyl(r, theta, phi)

spherical to cylindrical coordinate transform

# Arguments
- `r::Number` : spherical radius
- `theta::Number` : longitude, angle [0,2*pi]
- `phi::Number` : colatitude, angle [0.,pi]

# output
- `R, theta, z` : cylindrical coordiantes, R cylindrical radius, theta longitude, z, altitude

# Dependency
None
"""
function sph2cyl(r, theta, phi)
    R = r * sin(phi)
    z = r * cos(phi)
    return R, theta, z
end

    
#############################################################################
#### 3D INTEGRATION FUNCTIONS
#############################################################################
"""
    integ3D_spherical(fin,rmin,rmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)

3D integration of a complex function in spherical coordinates using cuhre

# Arguments
- `fin` : function to integrate fin(r,theta,phi,args)
- `rmin,rmax::Float64` : integration boundaries on r
- `args` : named tuple, container for the parameters of fin
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function integ3D_spherical(fin,rmin,rmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        r = rmin + x[1]*(rmax-rmin)
      	f[1],f[2] = reim(fin(r, x[2]*2.0*pi, x[3]*pi,args) * r^2* sin(x[3]*pi) *2.0* pi^2*(rmax-rmin))
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end
    

"""
    integ3D_galactic(fin,rmin,rmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)

3D integration ofa complex function in galactic coordinates (centered on the galaxy) using cuhre

# Arguments
- `fin` : function to integrate fin(r,ell,b,args)
- `rmin,rmax::Float64` : integration boundaries on r
- `args` : named tuple, container for the parameters of fin
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function integ3D_galactic(fin,rmin,rmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        r = rmin + x[1]*(rmax-rmin)
      	f[1],f[2] = reim(fin(r, x[2]*2.0*pi, x[3]*pi,args) * r^2* sin(pi/2.0-x[3]*pi) *2.0* pi^2*(rmax-rmin))
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end
    

"""
    integ3D_cylindrical(fin,rmin,rmax,zmin,zmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)

3D integration ofa complex function in galactic coordinates (centered on the galaxy) using cuhre

# Arguments
- `fin` : function to integrate fin(R,theta,z,args)
- `rmin,rmax::Float64` : integration boundaries on r
- `zmin,zmax::Float64` : integration boundaries on z
- `args` : named tuple, container for the parameters of fin
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function integ3D_cylindrical(fin,rmin,rmax,zmin,zmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        r = rmin + x[1]*(rmax-rmin)
        r = zmin + x[3]*(zmin-zmax)
      	f[1],f[2] = reim(fin(r, x[2]*2.0*pi, z,args) * r *2.0*pi*(rmax-rmin)*(zmax-zmin))
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end
    

"""
    integ3D_cartesian(fin,xmin,xmax,ymin,ymax,zmin,zmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)

3D integration of a complex function in cartesian coordinates using cuhre

# Arguments
- `fin` : function to integrate fin(x,y,z,args)
- `xmin,xmax::Float64` : integration boundaries on x
- `ymin,ymax::Float64` : integration boundaries on y
- `zmin,zmax::Float64` : integration boundaries on z
- `args` : named tuple, container for the parameters of fin
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function integ3D_cartesian(fin,xmin,xmax,ymin,ymax,zmin,zmax,args, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=100000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        xx = xmin + x[1]*(xmax-xmin)
        yy = ymin + x[2]*(ymax-ymin)
        zz = zmin + x[3]*(zmax-zmin)
      	f[1],f[2] = reim(fin(xx, yy, zz,args) * (xmax-xmin) * (ymax-ymin) * (zmax-zmin))
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end
    
    
#############################################################################
#### HERMITIAN PRODUCT FUNCTIONS
#############################################################################

"""
    dotprod_spherical(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in spherical coordinates using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(r,theta,phi,argsf)
- `gin` : function gin(r,theta,phi,argsg)
- `rmin, rmax::Float64` : minimum and maximum radius for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_spherical(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        r = rmin + x[1]*(rmax-rmin) 
      	f[1],f[2] = reim(fin(r, x[2]*2*pi, x[3]*pi,argsf) * conj(gin(r, x[2]*2*pi,x[3]*pi,argsg)) * r^2 * sin(x[3]*pi)* 2.0 *pi^2*(rmax-rmin))
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end


"""
    dotprod_galactic(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in galactic-like coordinates (centered on the galaxy) using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(r,ell,b,argsf)
- `gin` : function gin(r,ell,b,argsg)
- `rmin, rmax::Float64` : minimum and maximum radius for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_galactic(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        r = rmin + x[1]*(rmax-rmin)
        bb = x[3]*pi -pi/2.
      	f[1],f[2] = reim(fin(r, x[2]*2*pi, bb,argsf) * conj(gin(r, x[2]*2*pi, bb,argsg)) * r^2 * sin(pi/2.0 -x[3]*pi)* 2.0 *pi^2*(rmax-rmin))
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end

"""
    dotprod_cylindrical(fin, gin, Rmin::Float64, Rmax::Float64, zmin::Float64, zmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in cylindrical coordinates using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(R,theta,z,argsf)
- `gin` : function gin(R,theta,z,argsg)
- `Rmin, Rmax::Float64` : minimum and maximum radius for integration
- `zmin, zmax::Float64` : minimum and maximum z boundaries for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_cylindrical(fin, gin, Rmin::Float64, Rmax::Float64, zmin::Float64, zmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        R = Rmin + x[1]*(Rmax-Rmin)
        z = zmin + x[3]*(zmax-zmin)
      	f[1],f[2] = reim(fin(R, x[2]*2*pi, z,argsf) * conj(gin(R, x[2]*2*pi,z,argsg)) * R * (Rmax-Rmin) * (zmax-zmin) * 2 *pi)
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end

"""
    dotprod_cartesian(fin, gin, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, zmin::Float64, zmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in spherical coordinates using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(x,y,z,argsf)
- `gin` : function gin(x,y,z,argsg)
- `xmin, xmax::Float64` : boundaries along x for integration
- `ymin, ymax::Float64` : boundaries along y for integration
- `zmin, zmax::Float64` : boundaries along z for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_cartesian(fin, gin, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, zmin::Float64, zmax::Float64, argsf, argsg, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        xx = xmin + x[1] * (xmax-xmin)
        yy = ymin + x[2] * (ymin-ymax)
        zz = zmin + x[3] * (zmin-zmax)
      	f[1],f[2] = reim(fin(xx, yy, zz, argsf) * conj(gin(xx, yy, zz, argsg)) * (xmax-xmin) * (ymax-ymin) * (zmax-zmin))
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end



######## VECTORIAL FORM OF HERMITIAN PRODUCT

"""
    dotprod_spherical_vec(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in spherical coordinates using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(r,theta,phi,argsf)
- `gin` : function gin(r,theta,phi,argsg)
- `rmin, rmax::Float64` : minimum and maximum radius for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `nvec::Int=100` : number of vectors to use in cuhre
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_spherical_vec(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals=1000, maxevals=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        f[1,:] .= 1.0
        f[2,:] .= 1.0
        r = rmin .+ x[1,:]*(rmax-rmin) 
        for j = 1:size(x, 2) # loop over nvec
      	    f[1,j],f[2,j] = reim(fin(r[j], x[2,j]*2*pi, x[3,j]*pi,argsf) * conj(gin(r[j], x[2,j]*2*pi,x[3,j]*pi,argsg)) * r[j]^2 * sin(x[3,j]*pi)* 2.0 *pi^2*(rmax-rmin))
        end
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol,nvec=nvec, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end


"""
    dotprod_galactic_vec(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in galactic-like coordinates (centered on the galaxy) using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(r,ell,b,argsf)
- `gin` : function gin(r,ell,b,argsg)
- `rmin, rmax::Float64` : minimum and maximum radius for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `nvec::Int=100` : number of vectors to use in cuhre
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_galactic_vec(fin, gin, rmin::Float64, rmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals=1000, maxevals=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        f[1,:] .= 1.0
        f[2,:] .= 1.0
        r = rmin .+ x[1,:]*(rmax-rmin) 
        for j = 1:size(x, 2) # loop over nvec
      	    f[1,j],f[2,j] = reim(fin(r[j], x[2,j]*2*pi, x[3,j]*pi,argsf) * conj(gin(r[j], x[2,j]*2*pi,x[3,j]*pi,argsg)) * r[j]^2 * sin(pi/2.0-x[3,j]*pi)* 2.0 *pi^2*(rmax-rmin))
        end
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol,nvec=nvec, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end


"""
    dotprod_cylindrical_vec(fin, gin, Rmin::Float64, Rmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in cylindrical coordinates using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(R,theta,z,argsf)
- `gin` : function gin(R,theta,z,argsg)
- `Rmin, Rmax::Float64` : minimum and maximum radius for integration
- `zmin, zmax::Float64` : boundaries along z for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `nvec::Int=100` : number of vectors to use in cuhre
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_cylindrical_vec(fin, gin, Rmin::Float64, Rmax::Float64, zmin::Float64, zmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals=1000, maxevals=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        f[1,:] .= 1.0
        f[2,:] .= 1.0
        r = Rmin .+ x[1,:]*(Rmax-Rmin)
        zz = zmin .+ x[3,:]*(zmax-zmin)
        for j = 1:size(x, 2) # loop over nvec
      	    f[1,j],f[2,j] = reim(fin(r[j], x[2,j]*2*pi, zz[j],argsf) * conj(gin(r[j], x[2,j]*2*pi,zz[j],argsg)) * r[j] * 2.0 * pi * (Rmax-Rmin) * (zmin-zmax))
        end
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol,nvec=nvec, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end


"""
    dotprod_cartesian_vec(fin, gin, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, zmin::Float64, zmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals::Int=1000, maxevals::Int=1000000)

Computes the integral of fin . conj(gin) in cylindrical coordinates using Cuba cuhre algorithm's vectorial form for the integration

# Arguments
- `fin` : function fin(x,y,z,argsf)
- `gin` : function gin(x,y,z,argsg)
- `xmin, xmax::Float64` : boundaries along x for integration
- `ymin, ymax::Float64` : boundaries along y for integration
- `zmin, zmax::Float64` : boundaries along z for integration
- `argsf` : named tuple, container for the arguments of fin
- `argsg` : named tuple, container for the arguements of gin 
- `nvec::Int=100` : number of vectors to use in cuhre
- `atol::Float64=1e-12` : requested absolute accuracy for integration
- `rtol::Float64=1e-10` : requested relative accuracy for integration
- `minevals::Real=1000` : minimum number of evaluations for integration
- `maxevals::Real=1000000` : approximate maximum number of evaluations allowed

# Output
- `result` : result of the Hermitian product
- `error` : estimated error on the Hermitian product from cuhre
- `probability` : chi square probability that the error is not a reliable estimate of the true integration error (probability>0.95 indicates a potential problem)
- `neval` : number of integrand evaluations needed
- `fail` : error flag. 0 = desired accuracy was reached, -1=dimension out of range, 1 accuracy goal was not met.
- `nregions` : number of subregions needed

# Dependency
Cuba

# TODO
should offer possibility to integrate to infinity

"""
function dotprod_cartesian_vec(fin, gin, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, zmin::Float64, zmax::Float64, argsf, argsg, nvec::Int=100, atol::Float64=1e-12, rtol::Float64=1e-10, minevals=1000, maxevals=1000000)
    
    function integrand(x,f)
        # WARNING x is the NDarray of abscissa for integration in [0,1]^N
        f[1,:] .= 1.0
        f[2,:] .= 1.0
        xx = xmin .+ x[1,:]*(xmax-xmin) 
        yy = ymin .+ x[2,:]*(ymax-ymin) 
        zz = zmin .+ x[3,:]*(zmax-zmin) 
        for j = 1:size(x, 2) # loop over nvec
      	    f[1,j],f[2,j] = reim(fin(xx[j], yy[j], zz[j],argsf) * conj(gin(xx[j], yy[j], zz[j],argsg)) * (xmax-xmin) * (ymax-ymin) * (zmin-zmax))
        end
    end

    result, error, probability, neval, fail, nregions = cuhre(integrand, 3, 2, atol=atol, rtol=rtol,nvec=nvec, minevals=minevals,maxevals=maxevals)
    return result, error, probability, neval, fail, nregions
end

