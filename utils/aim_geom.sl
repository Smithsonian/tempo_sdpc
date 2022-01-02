#! /usr/bin/env slsh

require ("cmdopt");

% Geometry specified by Maxar and Intelsat:
variable Maxar_Geometry = struct
{
   lon_sc_deg = -91.0,    % -90.9 +/- 0.05 is Intelsat intended station
   bs_lat_deg =  33.71,   % boresight surface point
   bs_lon_deg = -91.65,
   % e.g. acos(bs_unit_vector) * 180/PI
   bs_angles_deg = [90.093, 95.447, 5.448]
};

private variable _DtoR   = PI/180.0;
private variable R_geo   = 42164.0;  % km
private variable R_earth = 6371.0;   % km (mean)

% WGS84 ellipsoid:
private variable R_a = 6378.137;  % km
private variable R_b = 6356.752;  % km
private variable _WGS84_ratio = R_a/R_b;

private define unit_vec (theta, phi)
{
   return [cos(phi)*sin(theta),
	   sin(phi)*sin(theta),
	   cos(theta)];
}

private define cross (a,b)
{
   return [ a[1] * b[2] - b[1] * a[2],
	  -(a[0] * b[2] - b[0] * a[2]),
	    a[0] * b[1] - b[0] * a[1]];
}

private define dot (a,b)
{
   return sum(a*b);
}

private define wgs84_xyz (lat_phi, lon_lam)
{
   variable a = R_a;
   variable b = R_b;
   variable h = 0.0;       % height above/below ellipsoid

   % In these expressions, phi = geodetic latitude, lambda = longitude
   variable phi = lat_phi * _DtoR;
   variable lam = lon_lam * _DtoR;
   variable cos_phi = cos(phi);
   variable sin_phi = sin(phi);

   variable N = a/hypot(cos_phi, (b/a) * sin_phi);
   variable x = (N + h) * cos_phi * cos(lam);
   variable y = (N + h) * cos_phi * sin(lam);
   variable z = (N * (b/a)^2 + h) * sin_phi;

   return [x,y,z];
}

private define geodetic_lat (geocentric_lat_deg)
{
   % tan(gencentric_lat) = (1-f)^2 tan(geodetic_lat)
   % where f = flattening = (a-b)/a
   % (1-f) = 1-(1-b/a) = b/a
   % so geodetic_lat = atan((a/b)^2 * tan(geocentric_lat))
   return atan(sqr(_WGS84_ratio) * tan(geocentric_lat_deg * _DtoR)) / _DtoR;
}

private define wgs84_xyz_geocentric (lat_deg, lon_deg)
{
   return wgs84_xyz (geodetic_lat (lat_deg), lon_deg);
}

private define sphere_xyz (lat_phi, lon_lam)
{
   % Unit vec wants theta, phi in the usual spherical coordinates
   variable   phi = lon_lam * _DtoR;
   variable theta = (90.0 - lat_phi) * _DtoR;

   return R_earth * unit_vec (theta, phi);
}

private define geocentric_lat_lon (r)
{
   variable r_plane = hypot (r[0], r[1]);
   variable lat_deg = atan2 (r[2], r_plane) / _DtoR;
   variable lon_deg = atan2 (r[1], r[0]) / _DtoR;
   return lat_deg, lon_deg;
}

private define coordinate_transform_matrices (lon_deg)
{
   variable sat_phi   = lon_deg * _DtoR;
   variable sat_theta = 0.5*PI;
   variable sat = R_geo * unit_vec(sat_theta, sat_phi);

   % spacecraft coordinate system unit vectors:
   % Geostationary spacecraft at longitude (lon_deg) above equator.
   % +Z toward the center of the earth
   % +Y toward south, parallel to earth rotation axis
   % +X along orbital velocity vector
   variable zhat = - unit_vec (sat_theta, sat_phi);
   variable yhat =   unit_vec (PI, 0);
   variable xhat =   cross (yhat, zhat);

   % coordinate transformation matrices between
   % WGS84 XYZ coordinates and satellite coordinate frame
   variable xyz_to_sat = _reshape ([xhat, yhat, zhat], [3,3]);
   variable sat_to_xyz = transpose(xyz_to_sat);

   return struct
     {
        sat_pos_xyz = sat,
        xyz_to_sat = xyz_to_sat,
        sat_to_xyz = sat_to_xyz
     };
}

private define compute_unit_vector_toward_point (lon_deg, pt)
{
   variable p = coordinate_transform_matrices (lon_deg);

   % vector in WGS84 (X,Y,Z) coordinates
   variable v = pt - p.sat_pos_xyz;
   variable v_u = v / hypot(v);

   % vector components in s/c system
   variable v_u_sc = p.xyz_to_sat # v_u;

   return v_u_sc;
}

private variable V = struct
{
   bs_u_xyz,     % unit vector along boresight in WGS84 XYZ coordinates
   sat_pos_xyz   % satellite position vector in WGS84 XYZ coordinates
};

% WGS84 XYZ position vector of the Earth surface point where
% the instrument boresight unit vector is aimed.
private define pos_vector_xyz (d0)
{
   return V.sat_pos_xyz + d0 * V.bs_u_xyz;
}

private define radial_distance_above_sphere (d0)
{
   return hypot(pos_vector_xyz (d0)) - R_earth;
}

private define radial_distance_above_ellipsoid (d0)
{
   variable pos = pos_vector_xyz (d0);
   variable lat_deg, lon_deg;
   (lat_deg, lon_deg) = geocentric_lat_lon (pos);
   variable ell = wgs84_xyz_geocentric (lat_deg, lon_deg);

   % distance along radius vector (not perpendicular to ellipsoid)
   return hypot(pos) - hypot(ell);
}

% Compute the (lon,lat) coordinates where the boresight unit vector
% from a geostationary satellite at the specified longitude
% intersects the earth's surface.
private define compute_pierce_point (lon_deg, bs_u_sc)
{
   variable is_sphere = qualifier_exists ("sphere");

   variable p = coordinate_transform_matrices (lon_deg);
   variable bs_u_xyz = p.sat_to_xyz # bs_u_sc;
   V.bs_u_xyz = bs_u_xyz;
   V.sat_pos_xyz = p.sat_pos_xyz;

   variable r;  % boresight vector
   variable d;  % scalar distance from satellite to boresight point

   % In S/C coordinates, the direction cosine for the z axis
   % is the cosine of the angle between the boresight direction
   % and the direction toward the earth's center (+z).
   % When the boresight intersects the earth, an upper limit
   % on the distance to the near-side intersection point is:
   variable dmax = R_geo * bs_u_sc[2];

   % The lower limit on the distance to the near-side
   % intersection point is the geostationary orbit height:
   variable dmin = R_geo - (is_sphere ? R_earth : R_a);

   variable radial_dist_func;
   if (is_sphere)
     radial_dist_func = &radial_distance_above_sphere;
   else
     radial_dist_func = &radial_distance_above_ellipsoid;

   variable num_loops = 0, converged = 0;
   variable d0 = dmin;
   variable d1 = dmax;
   variable h1, h0 = (@radial_dist_func)(d0);

   % Use secant method to find the root
   % (iterate linear interpolation to the zero)
   loop (64)
     {
        h1 = (@radial_dist_func)(d1);
        d = (d0 * h1 - d1 * h0) / (h1 - h0);
        if (abs (d1 - d) < 1.e-10*d)
          {
             converged = 1;
             break;
          }
        h0 = h1;
        d0 = d1;
        d1 = d;
        num_loops++;
     }

   if (converged == 0)
     throw ApplicationError, "*** iteration did not converge";

   %vmessage ("num_loops=$num_loops"$);

   % XYZ coordinates of earth surface intersection point
   r = pos_vector_xyz (d);

   % geocentric angular coordinates:
   variable lat, lon;
   (lat, lon) = geocentric_lat_lon (r);

   if (is_sphere == 0)
     {
        lat = geodetic_lat (lat);
     }

   return lon, lat;
}

private define compute_proj4_tpers_angles (bs_u_sc)
{
   % bs_u_sc has the boresight unit vector components
   % (e.g. direction cosines)
   % tpers azimuth zero is toward north pole (-y in S/C coordinate system)
   % S/C coordinates azimuth zero is along orbital velocity vector (+x)
   variable tilt = acos(bs_u_sc[2]) / _DtoR;
   variable az = 90.0 + atan2 (bs_u_sc[1], bs_u_sc[0]) / _DtoR;
   return tilt, az;
}

private define unit_vector_from_angles_deg (angles_deg)
{
   return cos(angles_deg * _DtoR);
}

private define angles_deg_from_unit_vector (bs_unit_vector)
{
   return acos(bs_unit_vector)/_DtoR;
}

private define maybe_convert_to_double (s)
{
   if (_typeof(s) == Double_Type)
     return s;
   return atof(s);
}

private define print_proj4_tpers_angles (bs_u)
{
   variable tilt, az;
   (tilt, az) = compute_proj4_tpers_angles (bs_u);
   vmessage ("Viewing angle parameters for proj4 tpers projection");
   vmessage ("%7.4f : [deg] tilt", tilt);
   vmessage ("%7.4f : [deg] azimuth", az);
}

private define print_boresight_components (bs_angles_deg)
{
   variable bs_u = unit_vector_from_angles_deg (bs_angles_deg);
   vmessage ("%12.4f, %12.4f, %12.4f  : direction angles [deg]", bs_angles_deg[0], bs_angles_deg[1], bs_angles_deg[2]);
   vmessage ("%12.4e, %12.4e, %12.4e  : direction cosines", bs_u[0], bs_u[1], bs_u[2]);
   print_proj4_tpers_angles (bs_u);
}

private define examine_geometry (geom)
{
   variable lon_sc_deg = geom.lon_sc_deg;
   variable bs_lat_deg = geom.bs_lat_deg;
   variable bs_lon_deg = geom.bs_lon_deg;
   variable bs_angles_deg = geom.bs_angles_deg;
   variable bs_u_calc, offset_rad, aim_lon, aim_lat, aim_pt;
   variable pt, bs_u;
   variable is_sphere = qualifier_exists ("sphere");

   if (lon_sc_deg == NULL)
     {
        vmessage ("*** Error: host spacecraft longitude not set");
        return 1;
     }

   vmessage ("===========================");
   vmessage ("Given:");
   vmessage ("%7.2f deg: spacecraft longitude", lon_sc_deg);
   if (bs_lat_deg != NULL && bs_lon_deg != NULL)
     {
        vmessage ("%7.2f deg: boresight aim point longitude", bs_lon_deg);
        vmessage ("%7.2f deg: boresight aim point latitude", bs_lat_deg);
     }
   if (bs_angles_deg != NULL)
     {
        vmessage ("boresight direction:");
        print_boresight_components (bs_angles_deg);
     }
   vmessage ("===========================");

   if (bs_lat_deg == NULL || bs_lon_deg == NULL)
     {
        bs_u = unit_vector_from_angles_deg (bs_angles_deg);
        (aim_lon, aim_lat) = compute_pierce_point (lon_sc_deg, bs_u ;; __qualifiers);
        vmessage ("computed aim point (lon, lat): = (%7.4f, %7.4f)", aim_lon, aim_lat);
     }
   else if (bs_angles_deg == NULL)
     {
        if (is_sphere)
          pt = sphere_xyz (bs_lat_deg, bs_lon_deg);
        else
          pt = wgs84_xyz (bs_lat_deg, bs_lon_deg);
        bs_u = compute_unit_vector_toward_point (lon_sc_deg, pt);
        bs_angles_deg = angles_deg_from_unit_vector (bs_u);
        vmessage ("boresight direction:");
        print_boresight_components (bs_angles_deg);
     }
   else
     {
        % This is the given boresight direction
        bs_u = unit_vector_from_angles_deg (bs_angles_deg);

        % This is the given boresight surface point
        if (is_sphere)
          pt = sphere_xyz (bs_lat_deg, bs_lon_deg);
        else
          pt = wgs84_xyz (bs_lat_deg, bs_lon_deg);

        % Check internal consistency:
        bs_u_calc = compute_unit_vector_toward_point (lon_sc_deg, pt);

        % For nearly collinear vectors, the cross product gives
        % a more robust measure of the angle between the vectors.
        offset_rad = asin(hypot(cross(bs_u, bs_u_calc)));
        vmessage ("Aim point and boresight vector are offset by :");
        vmessage ("%9.4f [urad]", offset_rad * 1.e6);
        vmessage ("%9.4f [km]", offset_rad * R_earth);

        % Computed surface point for the given boresight direction:
        (aim_lon, aim_lat) = compute_pierce_point (lon_sc_deg, bs_u;; __qualifiers);
        vmessage ("aim point for the given boresight direction:");
        vmessage ("lon: %9.4f deg (delta = %7.4f deg)", aim_lon, aim_lon - bs_lon_deg);
        vmessage ("lat: %9.4f deg (delta = %7.4f deg)", aim_lat, aim_lat - bs_lat_deg);

        % Computed boresight direction to given surface point:
        bs_angles_deg = angles_deg_from_unit_vector (bs_u_calc);
        vmessage ("boresight direction for given aim point:");
        print_boresight_components (bs_angles_deg);
     }

   return 0;
}

private define error_routine (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

private define usage()
{
   variable argv0 = __argv[0];
   variable s =
`
Usage: $argv0 [options]
       Options:
         -h|--help            Print this usage message
         -s|--lon_sc_deg LON  Spacecraft longitude [deg]
         -a|--bs_lat_deg LAT  Boresight aim point latitude [deg]
         -o|--bs_lon_deg LON  Boresight aim point longitude [deg]
         -b|--bs_angles VEC   Boresight direction angles [deg]
                              where VEC = "NULL" | "[alpha,beta,gamma]"
         -p|--bs_point        If present, compute boresight
                              aim point (lon,lat) coordinates
         -S|--sphere          Assume a spherical Earth

Default: WGS84 ellipsoid, geodetic latitude
`$;
   vmessage(s);
   exit(0);
}

define slsh_main()
{
   variable s = @Maxar_Geometry;
   variable bs_angles_str = NULL;
   variable compute_bs_point = 0;
   variable is_sphere = 0;

   variable c = cmdopt_new (&error_routine);
   c.add ("h|help", &usage);
   c.add ("s|lon_sc_deg", &s.lon_sc_deg; type="double");
   c.add ("a|bs_lat_deg", &s.bs_lat_deg; type="double");
   c.add ("o|bs_lon_deg", &s.bs_lon_deg; type="double");
   c.add ("b|bs_angles", &bs_angles_str; type="string");
   c.add ("p|bs_point", &compute_bs_point; inc);
   c.add ("S|sphere", &is_sphere; inc);
   variable __i = c.process (__argv, 1);

   if (_typeof(bs_angles_str) == String_Type)
     {
        s.bs_angles_deg = eval(bs_angles_str);
     }

   if (compute_bs_point)
     {
        s.bs_lat_deg = NULL;
        s.bs_lon_deg = NULL;
     }

   variable status;
   if (is_sphere)
     status = examine_geometry (s ; sphere);
   else
     status = examine_geometry (s);

   return status;
}
