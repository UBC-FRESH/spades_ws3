
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "spades_ws3",
  description = NA, #"insert module description here",
  keywords = NA, # c("insert key words here"),
  authors = c(person(c("First", "Middle"), "Last", email = "email@example.com", role = c("aut", "cre"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.5.9000", spades_ws3 = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "spades_ws3.Rmd"),
  reqdPkgs = list('R.utils'),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description
    defineParameter("basenames", "character", NA, NA, NA, 'MU baseneames to load'),
    defineParameter("base.year", 'numeric', 2015, NA, NA, 'base year of forest inventory data'),
    defineParameter("horizon", "numeric", 1, NA, NA, "ws3 simulation horizon (periods)"),
    defineParameter("tifPath", 'character', 'tif', NA, NA, desc = 'name of directory with tifs in inputs'),
    defineParameter("yearOfFirstHarvest", 'numeric', start(sim), NA, NA, "year to schedule first harvest"),
    defineParameter("scheduler.mode", "character", "optimize", NA, NA, "Switch between 'optimize' and 'areacontrol' harvest scheduler modes"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 10, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "landscape", objectClass = "RasterStack", desc = "stand age", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = 'landscape', objectClass = 'RasterStack', desc = 'raster stack of landscape attributes')
  )
))

## event types
#   - type `init` is required for initialization

doEvent.spades_ws3 = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "harvest")
      sim <- scheduleEvent(sim, start(sim), "spades_ws3", "grow")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "spades_ws3", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "spades_ws3", "save")
    },
    plot = {

      #plotFun(sim) # uncomment this, replace with object to plot
      # schedule future event(s)

      # e.g.,
      #sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "spades_ws3", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {

      sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "spades_ws3", "save")

      # schedule future event(s)

      # e.g.,
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "spades_ws3", "save")

      # ! ----- STOP EDITING ----- ! #
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
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {

  py$sys$path <- insert(py$sys$path, 1, file.path(modulePath(sim), currentModule(sim), "python"))
  py$sys$path <- insert(py$sys$path, 1, file.path(modulePath(sim), currentModule(sim), "python", "ws3"))
  py$basenames <- P(sim)$basenames
  #browser()
  py_run_file(file.path(modulePath(sim), currentModule(sim), "python", "spadesws3_params.py"))
  py$base_year <- P(sim)$base.year
  py$horizon <- P(sim)$horizon
  #browser()
  sim$fm <- py$bootstrap_forestmodel_kwargs()
  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot(sim$object)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}


updateAges <- function(sim) {
  #browser()
  year <- as.integer(time(sim) - start(sim) + P(sim, module = currentModule(sim))$base.year)
  files <- sapply(P(sim, module = currentModule(sim))$basenames,
                  function(bn) file.path(inputPath(sim),
                                         P(sim, module = currentModule(sim))$tifPath,
                                         bn,
                                         paste("inventory_", toString(year), ".tif", sep="")))
  rs.list <- sapply(files, stack) # one stack per MU
  rs.list <- rapply(rs.list,
                    function(rs, landscapeAge = sim$landscape$age) {
                      #browser()
                      rs[[2]] <- crop(sim$landscape$age, rs[[2]]) %>%
                        mask(., rs[[2]])
                      rs[[2]][is.nan(rs[[2]])] <- NA
                      return(rs)})
  mapply(writeRaster, rs.list, files, format='GTiff', overwrite=TRUE, datatype='INT4S')
  return(invisible(sim))
}

loadAges <- function(sim) {
  #browser()
  year <- as.integer(time(sim) - start(sim) + P(sim, module = currentModule(sim))$base.year)
  files <- sapply(P(sim, module = currentModule(sim))$basenames,
                  function(bn) file.path(inputPath(sim),
                                         P(sim, module = currentModule(sim))$tifPath,
                                         bn,
                                         paste("inventory_", toString(year), ".tif", sep="")))
  x <- sapply(files, raster, band=2)
   if (length(x) > 1) {
     x$fun <- mean
     x$na.rm <- TRUE
     browser()
     r <- do.call(mosaic, x)
     r[is.nan(r)] <- NA # replace NaN values with NA
   } else {
     r <- x[[1]]
   }
  names(x) <- NULL

  return(r)
}

applyHarvest <- function(sim) {
  #browser()
  year <- as.integer(time(sim) - start(sim) + P(sim, module = currentModule(sim))$base.year)
  py$base_year <- year
  sim$fm$base_year <- year
  updateAges(sim)
  masks <- paste(py$basenames, " 1 ? ?")
  areas <- c() # bogus placeholder (link this to user-defined input [list of values or function generating list of values])
  area.scale.factor <- 1. # bogus placeholder (will work)
  #browser()
  py$simulate_harvest(sim$fm,
                      list(py$basenames),
                      year,
                      P(sim, module=currentModule(sim))$scheduler.mode,
                      list(masks),
                      areas,
                      area.scale.factor) # run aspatial scheduler and allocate to pixels
  sim$landscape$age <- loadAges(sim)
  return(invisible(sim))
}

applyGrow <- function(sim) {
  sim$landscape$age <- sim$landscape$age + 1
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
