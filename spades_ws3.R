
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "spades_ws3",
  description = NA, #"insert module description here",
  keywords = NA, # c("insert key words here"),
  authors = c(person(c("Gregory", "Paradis"), "Last", email = "0@01101.io", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.5.9000", spades_ws3 = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "spades_ws3.Rmd"),
  reqdPkgs = list('R.utils', 'reticulate'),
  parameters = rbind(
    defineParameter("verbose", "numeric", 0, NA, NA, "console output verbosity level"),
    defineParameter("basenames", "character", NA, NA, NA, "MU baseneames to load"),
    defineParameter("horizon", "numeric", 1L, NA, NA, "ws3 simulation horizon (periods)"),
    defineParameter("base.year", "numeric", 2015L, NA, NA, "ws3 simulation base year"),
    defineParameter("scheduler.mode", "character", "optimize", NA, NA, "Switch between 'optimize' and 'areacontrol' harvest scheduler modes"),
    defineParameter("target.masks", "character", NULL, NA, NA, "Target masks (in '? ? ? ?' format). Only applicable if using 'areacontrol' scheduler mode."),
    defineParameter("target.areas", "numermic", NULL, NA, NA, "Target areas (ha).  Only applicable if using 'areacontrol' scheduler mode."),
    defineParameter("target.scalefactors", "numeric", NULL, NA, NA, "Target areas scale factors.  Only applicable if using 'areacontrol' scheduler mode."),
    defineParameter("mask.area.thresh", "numeric", 0., NA, NA, "Mask area threshold (for aggregation of bootstrapped masks).  Only applicable if using 'areacontrol' scheduler mode."),
    defineParameter("tifPath", 'character', 'tif', NA, NA, desc = 'name of directory with tifs in inputs'),
    defineParameter("yearOfFirstHarvest", 'numeric', start(sim), NA, NA, "year to schedule first harvest"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "landscape", objectClass = "RasterStack", desc = "stand age", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = 'landscape', objectClass = 'RasterStack', desc = 'raster stack of landscape attributes')
  )
))

## event types

doEvent.spades_ws3 = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "harvest")
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "grow")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "spades_ws3", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "spades_ws3", "save")
    },
    plot = {},
    save = {
      sim <- Save(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "spades_ws3", "save")
    },
    harvest = {
      sim <- applyHarvest(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "spades_ws3", "harvest")
    },
    grow = {
      sim <- applyGrow(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "spades_ws3", "grow")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

## event functions

Init <- function(sim) {
    #browser()
    library(R.utils)
    if (is.null(P(sim)$basenames)) stop(paste("'basenames' parameter value not specified in", currentModule(sim)))
    cmp <- grep(pattern = paste0(currentModule(sim), "$"), x = list.files(modulePath(sim))) %>%
           list.files(path = modulePath(sim), full.names = TRUE)[.] # current module path
    py$sys$path <- insert(py$sys$path, 1, file.path(cmp, "python"))
    py$sys$path <- insert(py$sys$path, 1, file.path(cmp, "python", "ws3"))
    py$basenames <- P(sim)$basenames
    py_run_file(file.path(cmp, "python", "spadesws3_params.py"))
    py$base_year <- P(sim)$base.year
    py$horizon <- P(sim)$horizon
    sim$fm <- py$bootstrap_forestmodel_kwargs()
    py$fm <- sim$fm
    return(invisible(sim))
}


Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}


plotFun <- function(sim) {
  return(invisible(sim))
}


updateAges <- function(sim, offset = 0) {
  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year)
  files1 <- sapply(P(sim)$basenames,
                   function(bn) file.path(inputPath(sim),
                                          P(sim)$tifPath,
                                          bn,
                                          paste("inventory_", toString(year), ".tif", sep="")))
  files2 <- sapply(P(sim)$basenames,
                   function(bn) file.path(inputPath(sim),
                                          P(sim)$tifPath,
                                          bn,
                                          paste("inventory_", toString(year+offset), ".tif", sep="")))
  #browser()
  rs.list <- sapply(files1, stack) # one stack per MU
  rs.list <- rapply(rs.list, 
                    function(rs) {
                    rs[[2]] <- crop(sim$landscape$age, rs[[2]]) %>% mask(., rs[[2]])
                    rs[[2]][is.nan(rs[[2]])] <- NA
                    return(rs)})
  mapply(writeRaster, rs.list, files2, format='GTiff', overwrite=TRUE, datatype='INT4S') 
  return(invisible(sim))
}


loadAges <- function(sim) {
  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year)
  files <- sapply(P(sim)$basenames,
                  function(bn) file.path(inputPath(sim),
                                         P(sim)$tifPath,
                                         bn,
                                         paste("inventory_", toString(year), ".tif", sep="")))
  x <- sapply(files, raster, band=2)
  if (length(x) > 1) {
    names(x)[1:2] <- c("x", "y") #from the raster pkg mosaic help. Needs x and y (!?)
    x$fun <- mean
    x$na.rm <- TRUE
    r <- do.call(mosaic, x)
    r[is.nan(r)] <- NA # replace NaN values with NA
  } else {
    r <- x[[1]]
  }
  names(x) <- NULL
  return(r)
}

applyHarvest <- function(sim) {
  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year)
  py$base_year <- year
  sim$fm$base_year <- year
  updateAges(sim)
  py$simulate_harvest(fm = sim$fm, 
                      basenames = P(sim)$basenames, 
                      year = year, 
                      mode = P(sim)$scheduler.mode, 
                      target_masks = P(sim)$target.masks, 
                      target_areas = P(sim)$target.areas,
                      target_scalefactors = P(sim)$target.scalefactors,
                      mask_area_thresh = P(sim)$mask.area.thresh,
                      verbose = P(sim)$verbose) 
  sim$landscape$age <- loadAges(sim)
  return(invisible(sim))
}


applyGrow <- function(sim) {
  sim$landscape$age <- sim$landscape$age + 1 
  updateAges(sim, offset=1)
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  return(invisible(sim))
}
