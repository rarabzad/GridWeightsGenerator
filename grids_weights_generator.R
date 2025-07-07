#' Generate grid-to-HRU intersection weights for regular or irregular grids
#'
#' @param ncfile   path to NetCDF file
#' @param hrufile  path to HRU shapefile (directory or .shp)
#' @param varnames character vector length-2 of NetCDF variable names (lon, lat), optional
#' @param dimnames character vector length-2 of NetCDF dimension names, optional
#' @param HRU_ID   name of the HRU ID field in the shapefile
#' @param plot     logical, whether to return a plotting function
#' @return         list(grid=sf, intersection=sf, centroids=sf, weights_txt=character, plot=function)
#' @import        sf ncdf4 geosphere dplyr sp
#' @export


grids_weights_generator <- function(ncfile, hrufile,
                                    varnames = NULL,
                                    dimnames = NULL,
                                    HRU_ID   = "HRU_ID",
                                    plot      = TRUE) {
  #-- 1) Dependencies and file checks -------------------------------------------
  required <- c('ncdf4','sf','geosphere','dplyr','sp','lwgeom','rmapshaper')
  miss     <- required[!sapply(required, requireNamespace, quietly=TRUE)]
  if (length(miss)) stop('Install missing packages: ', paste(miss, collapse=', '))
  if (!file.exists(ncfile))  stop('NetCDF not found: ', ncfile)
  if (!file.exists(hrufile)) stop('HRU file not found: ', hrufile)
  library(dplyr)
  sf::sf_use_s2(TRUE)
  
  #-- 2) Read NetCDF and build centroid matrices lonc, latc -----------------------
  nc <- tryCatch(ncdf4::nc_open(ncfile), error = function(e) stop('Failed to open NetCDF:', e$message))
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  if (!is.null(varnames))
  {
    if (length(varnames) != 2) stop("varnames must be a length-2 character vector")
    raw1 <- ncdf4::ncvar_get(nc, varnames[1])
    raw2 <- ncdf4::ncvar_get(nc, varnames[2])
    detect <- function(v) {
      if (grepl("lon", v, ignore.case = TRUE)) return("lon")
      if (grepl("lat", v, ignore.case = TRUE)) return("lat")
      u <- ncdf4::ncatt_get(nc, v, "units")$value
      if (!is.null(u)) {
        if (grepl("degree_east", u, ignore.case = TRUE)) return("lon")
        if (grepl("degree_north", u, ignore.case = TRUE)) return("lat")
      }
      NA_character_
    }
    kinds <- vapply(varnames, detect, character(1))
    if (any(is.na(kinds))) stop("Cannot auto-detect latitude vs. longitude variables")
    raw_lon <- if (kinds[1] == "lon") raw1 else raw2
    raw_lat <- if (kinds[1] == "lon") raw2 else raw1
    
    if (is.matrix(raw_lon) && is.matrix(raw_lat)) {
      if (!is.null(dimnames)) {
        perm_lon <- match(dimnames, sapply(nc$var[[varnames[kinds == "lon"]]]$dim, `[[`, "name"))
        perm_lat <- match(dimnames, sapply(nc$var[[varnames[kinds == "lat"]]]$dim, `[[`, "name"))
        if (any(is.na(c(perm_lon, perm_lat)))) stop("Provided dimnames don't match 2D var dims")
        lonc <- aperm(raw_lon, perm_lon)
        latc <- aperm(raw_lat, perm_lat)
      } else {
        lonc <- raw_lon
        latc <- raw_lat
      }
    } else if (is.vector(raw_lon) && is.vector(raw_lat)) {
      grid <- expand.grid(lon = raw_lon, lat = raw_lat)
      lonc <- matrix(grid$lon, nrow = length(raw_lat), ncol = length(raw_lon), byrow = FALSE)
      latc <- matrix(grid$lat, nrow = length(raw_lat), ncol = length(raw_lon), byrow = FALSE)
    } else {
      stop("varnames must point to either two 1D or two 2D variables")
    }
  } else {
    if (length(dimnames) != 2) stop("When varnames is NULL, dimnames must be length 2")
    dims <- nc$dim
    if (any(!dimnames %in% names(dims))) stop("Some dimnames not found in NetCDF dims")
    detect_dim <- function(dn)
    {
      var_dim <- dims[[dn]]
      if (grepl("lon", dn, ignore.case = TRUE)) return("lon")
      if (grepl("lat", dn, ignore.case = TRUE)) return("lat")
      if (!is.null(var_dim$units))
      {
        if (grepl("degree_east", var_dim$units, ignore.case = TRUE)) return("lon")
        if (grepl("degree_north", var_dim$units, ignore.case = TRUE)) return("lat")
      }
      NA_character_
    }
    kinds <- vapply(dimnames, detect_dim, character(1))
    if (any(is.na(kinds))) stop("Cannot detect which dimnames are lon/lat")
    lonv <- sort(dims[[dimnames[kinds == "lon"]]]$vals)
    latv <- sort(dims[[dimnames[kinds == "lat"]]]$vals)
    grid <- expand.grid(lon = lonv, lat = latv)
    lonc <- matrix(grid$lon, nrow = length(latv), ncol = length(lonv), byrow = FALSE)
    latc <- matrix(grid$lat, nrow = length(latv), ncol = length(lonv), byrow = FALSE)
  }
  
  if (!all(dim(lonc) == dim(latc))) stop("Centroid matrices lonc/latc dimensions mismatch")
  nr <- nrow(latc); ncg <- ncol(latc)
  
  #-- 3) Detect regular grid -----------------------------------------------------
  lon1d <- sort(unique(as.vector(lonc)))
  lat1d <- sort(unique(as.vector(latc)))
  is_regular <- (length(lon1d) > 1 && length(lat1d) > 1 &&
                   sd(diff(lon1d)) < .Machine$double.eps^0.5 &&
                   sd(diff(lat1d)) < .Machine$double.eps^0.5)
  
  #-- 4) Compute corner coordinates ----------------------------------------------
  if (is_regular)
  {
    dlon <- diff(lon1d)[1]; dlat <- diff(lat1d)[1]
    lon_edges <- c(lon1d[1] - dlon/2,
                   (lon1d[-length(lon1d)] + lon1d[-1]) / 2,
                   lon1d[length(lon1d)] + dlon/2)
    lat_edges <- c(lat1d[1] - dlat/2,
                   (lat1d[-length(lat1d)] + lat1d[-1]) / 2,
                   lat1d[length(lat1d)] + dlat/2)
    corner_lon <- matrix(lon_edges, nrow = nr+1, ncol = ncg+1, byrow = TRUE)
    corner_lat <- matrix(rev(lat_edges), nrow = nr+1, ncol = ncg+1, byrow = FALSE)
  } else {
    corner_lat <- matrix(NA_real_, nr+1, ncg+1)
    corner_lon <- corner_lat
    for (i in 1:(nr-1)) for (j in 1:(ncg-1)) {
      corner_lat[i+1, j+1] <- mean(latc[i:(i+1), j:(j+1)])
      corner_lon[i+1, j+1] <- mean(lonc[i:(i+1), j:(j+1)])
    }
    # Top/bottom edges
    for (j in 1:(ncg-1)) {
      top_mid <- c(mean(lonc[1, j:(j+1)]), mean(latc[1, j:(j+1)]))
      below_mid <- c(mean(lonc[2, j:(j+1)]), mean(latc[2, j:(j+1)]))
      br_top <- geosphere::bearing(below_mid, top_mid); d_top <- geosphere::distGeo(below_mid, top_mid)/2
      ext_top <- geosphere::destPoint(top_mid, br_top, d_top)
      corner_lon[1, j+1] <- ext_top[1]; corner_lat[1, j+1] <- ext_top[2]
      
      bot_mid <- c(mean(lonc[nr, j:(j+1)]), mean(latc[nr, j:(j+1)]))
      above_mid <- c(mean(lonc[nr-1, j:(j+1)]), mean(latc[nr-1, j:(j+1)]))
      br_bot <- geosphere::bearing(above_mid, bot_mid); d_bot <- geosphere::distGeo(above_mid, bot_mid)/2
      ext_bot <- geosphere::destPoint(bot_mid, br_bot, d_bot)
      corner_lon[nr+1, j+1] <- ext_bot[1]; corner_lat[nr+1, j+1] <- ext_bot[2]
    }
    # Left/right edges
    for (i in 1:(nr-1)) {
      left_mid <- c(mean(lonc[i:(i+1), 1]), mean(latc[i:(i+1), 1]))
      inl      <- c(mean(lonc[i:(i+1), 2]), mean(latc[i:(i+1), 2]))
      br_l     <- geosphere::bearing(inl, left_mid); d_l <- geosphere::distGeo(inl, left_mid)/2
      ext_l    <- geosphere::destPoint(left_mid, br_l, d_l)
      corner_lon[i+1, 1] <- ext_l[1]; corner_lat[i+1, 1] <- ext_l[2]
      
      right_mid <- c(mean(lonc[i:(i+1), ncg]), mean(latc[i:(i+1), ncg]))
      inr       <- c(mean(lonc[i:(i+1), ncg-1]), mean(latc[i:(i+1), ncg-1]))
      br_r      <- geosphere::bearing(inr, right_mid); d_r <- geosphere::distGeo(inr, right_mid)/2
      ext_r     <- geosphere::destPoint(right_mid, br_r, d_r)
      corner_lon[i+1, ncg+1] <- ext_r[1]; corner_lat[i+1, ncg+1] <- ext_r[2]
    }
    # Outer corners
    corner_list <- list(
      NW=list(pt=c(lonc[1,1],latc[1,1]), m=c(mean(c(lonc[2,1],lonc[1,2])), mean(c(latc[2,1],latc[1,2]))), idx=c(1,1)),
      NE=list(pt=c(lonc[1,ncg],latc[1,ncg]), m=c(mean(c(lonc[2,ncg],lonc[1,ncg-1])), mean(c(latc[2,ncg],latc[1,ncg-1]))), idx=c(1,ncg+1)),
      SW=list(pt=c(lonc[nr,1],latc[nr,1]), m=c(mean(c(lonc[nr-1,1],lonc[nr,2])), mean(c(latc[nr-1,1],latc[nr,2]))), idx=c(nr+1,1)),
      SE=list(pt=c(lonc[nr,ncg],latc[nr,ncg]), m=c(mean(c(lonc[nr-1,ncg],lonc[nr,ncg-1])), mean(c(latc[nr-1,ncg],latc[nr,ncg-1]))), idx=c(nr+1,ncg+1))
    )
    for (info in corner_list) {
      ext_pt <- geosphere::destPoint(info$pt, geosphere::bearing(info$m, info$pt), geosphere::distGeo(info$m, info$pt))
      corner_lon[info$idx[1], info$idx[2]] <- ext_pt[1]
      corner_lat[info$idx[1], info$idx[2]] <- ext_pt[2]
    }
  }
  
  #-- 5) Build grid-cell polygons ------------------------------------------------
  total <- nr * ncg
  polys <- vector("list", total)
  ids   <- rows <- cols <- integer(total)
  idx <- 1
  for (i in 1:nr) for (j in 1:ncg)
  {
    coords <- matrix(c(
      corner_lon[i, j], corner_lat[i, j],
      corner_lon[i, j+1], corner_lat[i, j+1],
      corner_lon[i+1, j+1], corner_lat[i+1, j+1],
      corner_lon[i+1, j], corner_lat[i+1, j],
      corner_lon[i, j], corner_lat[i, j]
    ), ncol=2, byrow=TRUE)
    polys[[idx]] <- sf::st_polygon(list(coords))
    ids[idx]     <- (j - 1) * nr + i-1
    rows[idx]    <- i; cols[idx] <- j; idx <- idx + 1
  }
  cent_sf <- sf::st_as_sf(
    data.frame(x = as.vector(lonc), y = as.vector(latc)),
    coords = c('x','y'), crs = 4326
  )
  grid_sf <- sf::st_sf(
    Cell_ID  = ids,       # 1-based ID
    geometry = sf::st_sfc(polys),
    crs      = 4326
  )
  grid_sf<-grid_sf[unlist(st_contains(grid_sf,cent_sf)),]
  grid_sf$Cell_ID<-0:(nrow(grid_sf)-1)
  cent_sf$Cell_ID<-0:(nrow(cent_sf)-1)
  
  
  #-- 6) Intersection with HRUs & weight calculation -----------------------------
  sf::sf_use_s2(FALSE)
  hru_sf <- sf::st_read(hrufile, quiet=TRUE)[, HRU_ID] %>% sf::st_transform(sf::st_crs(grid_sf))
  inter_sf <- sf::st_intersection(grid_sf, hru_sf %>% mutate(HRU_ID = row_number())) %>%
    group_by(HRU_ID) %>%
    mutate(weight = as.numeric(sf::st_area(geometry) / sum(sf::st_area(geometry), na.rm=TRUE))) %>%
    ungroup()
  
  #-- 7) Prepare outputs ---------------------------------------------------------
  wt <- inter_sf %>% sf::st_set_geometry(NULL) %>% dplyr::select(HRU_ID, Cell_ID, weight)
  weights_txt <- c(
    ':GridWeights',
    paste0(':NumberHRUs       ', n_distinct(wt$HRU_ID)),
    paste0(':NumberGridCells ', nr * ncg),
    '# HRU_ID\tCell_ID\tweight',
    apply(wt, 1, paste, collapse='\t'),
    ':EndGridWeights'
  )
  
  #-- Optional plot function ----------------------------------------------------
  plot_fn <- NULL
  if (plot) {
    plot_fn <- function(filename = "map.pdf") {
      pdf(filename, width = 8, height = 8)
      plot(as_Spatial(grid_sf), col = 'lightgrey', main = "Overlay Map")
      plot(as_Spatial(cent_sf), add = TRUE, col = 'red', pch = 19, cex = 0.5)
      plot(as_Spatial(rmapshaper::ms_simplify(hru_sf, keep = 0.05, keep_shapes = TRUE)),
           add = TRUE, border = 'blue')
      dev.off()
      normalizePath(filename)
    }
  }
  
  return(list(
    grid_sf     = grid_sf,
    hru_sf      = inter_sf,
    centroids   = cent_sf,
    weights_txt = weights_txt,
    plot_path   = plot_fn
  ))
}
