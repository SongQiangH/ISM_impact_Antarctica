undef("cal_rws")
; The function for calculating RWS (Rossby Wave Source) is sourced from: https://github.com/luruia/Tropic_CAA/blob/main/figs6%2B7.ncl
; Zhu, Z., Lu, R., Yu, B. et al. A moderator of tropical impacts on climate in Canadian Arctic Archipelago during boreal summer. Nat Commun 15, 8644 (2024). https://doi.org/10.1038/s41467-024-53056-0
function cal_rws(u_all, v_all)     
    local dims,lat,lon,nlat,mlon,uvmsg,dv,vr,uvd,uvdu,uvdv,dv_d,vr_d,dv_clm,vr_clm,uvd_clm, \
            uvdu_clm,uvdv_clm,dv_d_clm,vr_d_clm,f,f_2d,f_3d,abo_vor,abo_vor_clm,dlon,dlat, \
            dvordx,dvordy,dvordx_clm,dvordy_clm,nl,ml,dx,dy,S1,S2,S3,S4,RWS,u_clim,v_clim
    begin
        u_clim = dim_avg_n_Wrap(u_all, 0)
        v_clim = dim_avg_n_Wrap(v_all, 0)
        u_anom = dim_rmvmean_n_Wrap(u_all, 0)
        v_anom = dim_rmvmean_n_Wrap(v_all, 0)
        dims = dimsizes(u_anom)
        uvmsg = u_anom@_FillValue
        lat = u_anom&latitude
        lon = u_anom&longitude
        nlat = dimsizes(lat)
        mlon = dimsizes(lon)

        dv_anom = uv2dv_cfd(u_anom, v_anom, lat, lon, 3)
        vr_anom = uv2vr_cfd(u_anom, v_anom, lat, lon, 3)
        copy_VarMeta(u_anom, dv_anom)
        copy_VarMeta(u_anom, vr_anom)

        uvd_anom = dv2uvF_Wrap(dv_anom)   
        uvdu_anom = uvd_anom(0,:,:,:)
        uvdv_anom = uvd_anom(1,:,:,:)
        copy_VarMeta(u_anom, uvdu_anom)
        copy_VarMeta(u_anom, uvdv_anom)

        dv_d_anom = uv2dv_cfd(uvdu_anom, uvdv_anom, lat, lon, 1)
        vr_d_anom = uv2vr_cfd(uvdu_anom, uvdv_anom, lat, lon, 1)
        copy_VarMeta(u_anom, dv_d_anom)
        copy_VarMeta(u_anom, vr_d_anom)


        dv_clm = uv2dv_cfd(u_clim, v_clim, lat, lon, 3)
        vr_clm = uv2vr_cfd(u_clim, v_clim, lat, lon, 3)
        copy_VarMeta(u_clim, dv_clm)
        copy_VarMeta(u_clim, vr_clm)


        uvd_clm = dv2uvF_Wrap(dv_clm)   
        uvdu_clm = uvd_clm(0,:,:)
        uvdv_clm = uvd_clm(1,:,:)
        copy_VarMeta(u_clim, uvdu_clm)
        copy_VarMeta(u_clim, uvdv_clm)

        dv_d_clm = uv2dv_cfd(uvdu_clm, uvdv_clm, lat, lon, 1)
        vr_d_clm = uv2vr_cfd(uvdu_clm, uvdv_clm, lat, lon, 1)
        copy_VarMeta(u_clim, dv_d_clm)
        copy_VarMeta(u_clim, vr_d_clm)

        f = coriolis_param(lat)
        f_2d = conform(u_clim, f, (/0/))
        f_3d = conform(u_anom, f, (/1/))
        abo_vor_anom = vr_anom + f_3d   
        abo_vor_clm = vr_clm + f_2d   

        ;>---------------------------------------------------------------------------------------<
        ;>---------------------------------------------------------------------------------------<

        dlon = (lon(2)-lon(1)) * 0.0174533
        dlat = (lat(2)-lat(1)) * 0.0174533

        dvordx_anom = new ( (/dims(0), nlat, mlon/), "float", uvmsg )
        dvordy_anom = new ( (/dims(0), nlat, mlon/), "float", uvmsg )
        dvordx_clm = new ( (/nlat, mlon/), "float", uvmsg )
        dvordy_clm = new ( (/nlat, mlon/), "float", uvmsg )

        gradsf(abo_vor_clm,dvordx_clm,dvordy_clm)
        gradsf(abo_vor_anom,dvordx_anom,dvordy_anom)


        ;>---------------------------------------------------------------------------------------<
        ;>---------------------------------------------------------------------------------------<

        uvdu_clm := conform(uvdu_anom, uvdu_clm, (/1, 2/))
        uvdv_clm := conform(uvdv_anom, uvdv_clm, (/1, 2/))
        abo_vor_clm := conform(abo_vor_anom, abo_vor_clm, (/1, 2/))
        dvordx_clm := conform(dvordx_anom, dvordx_clm, (/1, 2/))
        dvordy_clm := conform(dvordy_anom, dvordy_clm, (/1, 2/))
        dv_d_clm := conform(dvordy_anom, dv_d_clm, (/1, 2/))
        
        S1 = - uvdu_anom * dvordx_clm - uvdv_anom * dvordy_clm
        S2 = - abo_vor_clm * dv_d_anom 
        S3 = - uvdu_clm * dvordx_anom - uvdv_clm * dvordy_anom
        S4 = - abo_vor_anom * dv_d_clm

        S1 = S1 * 10^11 ;−v′χ·∇ζ̄ 
        S2 = S2 * 10^11 ;−ζ̄ ∇·χ′ 
        S3 = S3 * 10^11 ;−v̄χ·∇ζ′
        S4 = S4 * 10^11 ;−ζ′∇·v̄χ

        RWS = S1 + S2 + S3 + S4

        copy_VarCoords(u_anom, S1)
        copy_VarCoords(u_anom, S2)
        copy_VarCoords(u_anom, S3)
        copy_VarCoords(u_anom, S4)
        copy_VarCoords(u_anom, RWS)

        return [/RWS, S1, S2, S3, S4/]
    end

fu = addfile("/public/home/songqh/project/ISM_impact_Antarctic/data/uwnd_jja_200hPa.nc", "r")
fv = addfile("/public/home/songqh/project/ISM_impact_Antarctic/data/vwnd_jja_200hPa.nc", "r")
u_all = fu->u  
v_all = fv->v
lat = fu->latitude
lon = fu->longitude
cache = cal_rws(u_all,v_all)
rws = cache[0]
s1 = cache[1]
s2 = cache[2]
s3 = cache[3]
s4 = cache[4]

system("rm -rf /public/home/songqh/project/ISM_impact_Antarctic/data/rws_output.nc")
fout = addfile("/public/home/songqh/project/ISM_impact_Antarctic/data/rws_output.nc", "c")

year = ispan(1979, 2024, 1)
year!0 = "year"
fout->year = year
fout->lat = lat
fout->lon = lon

rws!0 = "year"
rws!1 = "lat"
rws!2 = "lon"
rws@long_name = "Rossby Wave Source"
rws@units = "10^-11 s-2" 
rws@_FillValue = rws@_FillValue  
fout->rws = rws

s1!0 = "year"
s1!1 = "lat"
s1!2 = "lon"
s1@long_name = "RWS Term S1: −v′χ·∇ζ̄ "
s1@units = "10^-11 s-2" 
s1@_FillValue = s1@_FillValue  
fout->s1 = s1

s2!0 = "year"
s2!1 = "lat"
s2!2 = "lon"
s2@long_name = "RWS Term S2: −ζ̄ ∇·χ′ "
s2@units = "10^-11 s-2" 
s2@_FillValue = s2@_FillValue  
fout->s2 = s2

s3!0 = "year"
s3!1 = "lat"
s3!2 = "lon"
s3@long_name = "RWS Term S3: −v̄χ·∇ζ′"
s3@units = "10^-11 s-2" 
s3@_FillValue = s3@_FillValue  
fout->s3 = s3

s4!0 = "year"
s4!1 = "lat"
s4!2 = "lon"
s4@long_name = "RWS Term S4: −ζ′∇·v̄χ"
s4@units = "10^-11 s-2" 
s4@_FillValue = s4@_FillValue  
fout->s4 = s4