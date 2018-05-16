

ninoplot <-
  function (object,
            axes = c(1, 2),
            geom = "point",
            layers = c("species", "sites", "biplot", "centroids"),
            arrows = TRUE,
            legend.position = "right",
            palettage = FALSE,
            title = NULL,
            subtitle = NULL,
            caption = NULL,
            ylab = NA,
            xlab = NA,
            const,
            mypal = "PiYG",
            quali_sup = NULL,
            ...)

#Settings used to dev the function

# axes = c(1, 2)
# object <- rda_trait
# geom = "point"
# layers = c("species", "sites", "biplot", "centroids")
# arrows = TRUE
# legend.position = "right"
# title = NULL
# subtitle = NULL
# caption = NULL
# ylab
# xlab
# const
  {
    library(ggvegan)
    axes <- rep(axes, length.out = 2L)
    axnames <- paste0("RDA", axes)
    exp_var <- summary(object)$concont$importance[2, axes]
    axnames <- paste0(axnames, " (", round(100 * exp_var, 2), "%)")
    # axnames check
    obj <- fortify(object, axes = axes)
    LAYERS <- levels(obj$Score)
    vars <- ggvegan:::getDimensionNames(obj)
    # geom <- match.arg(geom)
    point <- TRUE
    if (isTRUE(all.equal(geom, "text"))) {
      point <- FALSE
    }
    obj <- obj[obj$Score %in% layers, , drop = FALSE]
    plt <- ggplot()
    if (isTRUE(arrows)) {
      take <- "sites"
    }
    if (!isTRUE(arrows)) {
      take <- c("species", "sites")
    }
    want <- obj$Score %in% take
    mypoints <- obj[want, , drop = FALSE] #hack here
    if(!is.null(quali_sup)){ #new argument to add thematic coloration
      mypoints$Morphotype <- quali_sup #hack here
    }
    
    
    if (point) {
      if(palettage == TRUE){
        plt <- plt + geom_point(data = mypoints,
                                # hack here
                                aes_string(
                                  x = vars[1],
                                  y = vars[2],
                                  # shape = "Score",
                                  colour = "Morphotype"
                                )) +
          scale_color_brewer(palette = mypal)
      }
      else {
        plt <- plt + geom_point(data = mypoints,
                                # hack here
                                aes_string(
                                  x = vars[1],
                                  y = vars[2],
                                  # shape = "Score",
                                  colour = "Morphotype"
                                )) #hack here
      }

      #HERE
    }
    # plt
    #
    # str(mypoints)
    # obj[want, , drop = FALSE] %>% ncol
    # obj[want, , drop = FALSE] %>% nrow
    # data_standar$morphotype %>% length
    #
    if (!point) {
      plt <- plt + geom_text(
        data = mypoints,
        aes_string(
          x = vars[1],
          y = vars[2],
          label = "Label",
          colour = "Morphotype"
        ),
        size = 3
      )
    }
    # plt
    if (isTRUE(arrows)) {
      want <- obj$Score == "species"
      pdat <- obj[want, , drop = FALSE]
      col <- "black"
      want <- obj[["Score"]] == "species"
      plt <- plt + geom_segment(
        data = pdat,
        aes_string(
          x = 0,
          y = 0,
          xend = vars[1],
          yend = vars[2]
        ),
        arrow = arrow(length = unit(0.2,
                                    "cm")),
        colour = col
      )
      pdat[, vars] <- 1.1 * pdat[, vars, drop = FALSE]
      plt <- plt + geom_text(data = pdat,
                             aes_string(x = vars[1],
                                        y = vars[2], label = "Label"),
                             size = 4)
    }
    if (all(c("biplot", "centroids") %in% LAYERS)) {
      want <- obj$Score == "biplot"
      tmp <- obj[want,]
      obj <- obj[!want,]
      bnam <- tmp[, "Label"]
      cnam <- obj[obj$Score == "centroids", "Label"]
      obj <- rbind(obj, tmp[!bnam %in% cnam, , drop = FALSE])
    }
    if (any(want <- obj$Score == "constraints")) {
      if (point) {
        plt <- plt + geom_point(data = obj[want, , drop = FALSE],
                                aes_string(x = vars[1], y = vars[2]))
      }
      else {
        plt <- plt + geom_text(data = obj[want, , drop = FALSE],
                               aes_string(
                                 x = vars[1],
                                 y = vars[2],
                                 label = "Label"
                               ))
      }
    }
    if (any(want <- obj$Score == "biplot")) {
      if (length(layers) > 1) {
        mul <- ggvegan:::arrowMul(obj[want, vars, drop = FALSE], obj[!want,
                                                           vars, drop = FALSE])
        obj[want, vars] <- mul * obj[want, vars]
      }
      col <- "red"
      plt <- plt + geom_segment(
        data = obj[want, , drop = FALSE],
        aes_string(
          x = 0,
          y = 0,
          xend = vars[1],
          yend = vars[2]
        ),
        arrow = arrow(length = unit(0.2, "cm")),
        colour = col
      )
      obj[want, vars] <- 1.1 * obj[want, vars]
      plt <- plt + geom_text(data = obj[want, , drop = FALSE],
                             aes_string(x = vars[1], y = vars[2], label = "Label"), colour = "red")
    }
    if (any(want <- obj$Score == "centroids")) {
      plt <- plt + geom_text(
        data = obj[want, , drop = FALSE],
        aes_string(x = vars[1], y = vars[2], label = "Label"),
        colour = "navy"
      )
    }
    if (is.na(xlab)) {
      xlab <- axnames[1]
    }
    if (is.na(ylab)) {
      ylab <- axnames[2]
    }
    plt <-
      plt + labs(
        x = xlab,
        y = ylab,
        title = title,
        subtitle = subtitle,
        caption = caption
      )
    plt <- plt + coord_fixed(ratio = 1)
    plt <- plt + theme(legend.position = legend.position)
    plt
  }
