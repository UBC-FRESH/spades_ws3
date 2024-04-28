
# Everything in this file gets sourced during simInit, and all functions 
# and objects are put into the simList. To use objects, use sim$xxx, and 
# are thus globally available to all modules. Functions can be used without 
# sim$ as they are namespaced, like functions in R packages. If exact location
# is required, functions will be: sim$<moduleName>$FunctionName.
library(dplyr)

defineModule(sim, list(
  name = "spades_ws3",
  description = "", 
  keywords = c(""), 
  authors = c(person(c("Gregory"), "Paradis", email = "gregory.paradis@ubc.ca", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.0.0", spades_ws3 = "2.0.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "spades_ws3.Rmd"),
  reqdPkgs = list("R.utils", "reticulate"),
  parameters = rbind(
    defineParameter("verbose", "numeric", 0, NA, NA, "console output verbosity level"),
    defineParameter("horizon", "numeric", 1L, NA, NA, "ws3 simulation horizon (periods)"),
    defineParameter("period.length", "numeric", 1L, NA, NA, "ws3 simulation period length"),
    defineParameter("base.year", "numeric", 2015L, NA, NA, "ws3 simulation base year"),
    defineParameter("model.name", "character", NA, NA, NA, "ws3 model name"),
    defineParameter("model.path", "character", NA, NA, NA, "ws3 model path"),
    defineParameter("shp.path", "character", NA, NA, NA, "ws3 initial inventory Shapefile file path"),
    defineParameter("tif.path", "character", NA, NA, NA, "ws3 GeoTIFF file path"),
    defineParameter("theme.cols", "character", NA, NA, NA, "list of theme column names"),
    defineParameter("scheduler.mode", "character", "optimize", NA, NA, "Switch between 'optimize' and 'areacontrol' harvest scheduler modes"),
    #defineParameter("target.masks", "character", NULL, NA, NA, "Target masks (in '? ? ? ?' format). Only applicable if using 'areacontrol' scheduler mode."),
    #defineParameter("target.areas", "numermic", NULL, NA, NA, "Target areas (ha).  Only applicable if using 'areacontrol' scheduler mode."),
    #defineParameter("target.scalefactors", "numeric", NULL, NA, NA, "Target areas scale factors.  Only applicable if using 'areacontrol' scheduler mode."),
    #defineParameter("mask.area.thresh", "numeric", 0., NA, NA, "Mask area threshold (for aggregation of bootstrapped masks).  Only applicable if using 'areacontrol' scheduler mode."),
    #defineParameter("tifPath", 'character', 'tif', NA, NA, desc = 'name of directory with tifs in inputs'),
    #defineParameter("yearOfFirstHarvest", 'numeric', start(sim), NA, NA, "year to schedule first harvest"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "landscape", objectClass = "RasterStack", desc = "raster stack of landscape attributes (age, hashed development type key)", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "landscape", objectClass = "RasterStack", desc = "raster stack of landscape attributes (age, hashed development type key)")
  )
))

## event types

doEvent.spades_ws3 = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "fire")
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "actions")
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "grow")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "spades_ws3", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "spades_ws3", "save")
    },
    plot = {},
    save = {
      sim <- Save(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "spades_ws3", "save")
    },
    fire = {
      sim <- applyFire(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "spades_ws3", "fire")
    },
    actions = {
      sim <- applyActions(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "spades_ws3", "actions")
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
  library(raster)
  # cmp <- grep(pattern = paste0(currentModule(sim), "$"), x = list.files(modulePath(sim))) %>%
  #             list.files(path = modulePath(sim), full.names = TRUE)[.] # current module path
  #py_run_file(file.path(cmp, "spadesws3.py"))
  #py_run_file("spadesws3.py")
  sim$spadesws3 <- import_from_path("spadesws3")
  sim$spadesws3$init_action_decode(P(sim)$action.decode)
  #lapply(P(sim)$action_decode, sim$spadesws3$add_raster_action)
  #tif.filename <- paste("inventory_", toString(P(sim)$base.year), ".tif", sep="")
  sim$spadesws3$import_inventory(shp_path=P(sim)$shp.path, 
                                 tif_path=P(sim)$tif.path, 
                                 tif_filename="inventory_init.tif",
                                 theme_cols=P(sim)$theme.cols)
  # sim$hdtk.decode <- reticulate::dict(py$import_inventory(shp_path=P(sim)$shp.path, 
  #                                                         tif_path=P(sim)$tif.path, 
  #                                                         tif_filename="inventory_init.tif",
  #                                                         theme_cols=P(sim)$theme.cols),
  #                                     convert=TRUE)
  sim$landscape <- stack(file.path(P(sim)$tif.path, "inventory_init.tif"))
  names(sim$landscape) <- c('hdtk', 'age', 'blockid', 'actions')
  return(invisible(sim))
}

Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}


plotFun <- function(sim) {
  return(invisible(sim))
}

applyFire <- function(sim) {
  print('sum area before fire', str(cellStats(sim$landscape$age, sum)))
  # BOGUS fire (simulates other modules requesting ws3 actions)
  sim$landscape$actions <- NAvalue(sim$landscape$actions) # reset
  maxage <- 150 - (10 * time(sim)) # burn oldest pixels first, in 10 year chunks
  print('maxage', str(maxage))
  sim$landscape$actions[sim$landscape$age > maxage] <- 1 # 1 is 'fire' set in P(sim)$action.decode
  print('sum area after fire', str(cellStats(sim$landscape$age, sum)))
  return(invisible(sim))
}

applyActions <- function(sim) {
  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year)
  tif.filename <- paste("inventory_", toString(year), ".tif", sep="")
  sim$foo1 <- c(cellStats(sim$landscape$age, sum)) # debug
  writeRaster(sim$landscape, file.path(P(sim)$tif.path, tif.filename), 
              format='GTiff', overwrite = TRUE, datatype = "INT4S")
  sim$spadesws3$bootstrap_forestmodel(model_name=P(sim)$model.name,
                                      model_path=P(sim)$model.path,
                                      base_year=year,
                                      tif_path=P(sim)$tif.path,
                                      tif_filename=tif.filename)
  sim$spadesws3$simulate_actions(scheduler_mode=P(sim)$scheduler.mode,
                                 targets=P(sim)$scheduler.targets)
  #sim$spadesws3$forestmodel$forestraster$allocate_schedule()
  sim$landscape <- stack(file.path(P(sim)$tif.path, tif.filename))
  names(sim$landscape) <- c('hdtk', 'age', 'blockid', 'actions')
  sim$foo2 <- c(cellStats(sim$landscape$age, sum)) # debug
  return(invisible(sim))
}


applyGrow <- function(sim) {
  print('sum area before grow', str(cellStats(sim$landscape$age, sum)))
  sim$foo3 <- c(cellStats(sim$landscape$age, sum))
  sim$landscape$age <- sim$landscape$age + 1 
  year <- as.integer(time(sim) - start(sim) + P(sim)$base.year + 1)
  tif.filename <- paste("inventory_", toString(year), ".tif", sep="")
  writeRaster(sim$landscape, file.path(P(sim)$tif.path, tif.filename), 
              format='GTiff', overwrite = TRUE, datatype = "INT4S")
  print('sum area after grow', str(cellStats(sim$landscape$age, sum)))
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  return(invisible(sim))
}
