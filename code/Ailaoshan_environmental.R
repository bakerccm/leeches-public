# R code to process environmental data

library("here")
library("tidyverse")
library("readxl")

# get data

    env.data <- read_excel(here("data", "environmental_variables_20180801.xlsx"), sheet = "variables", col_names = TRUE)

# rename polygon id to be consistent with other existing data

    env.data <- env.data %>%
        dplyr::rename(Polygon_ID = polygon_ID) %>%
        mutate(Polygon_ID = as.character(Polygon_ID))

# where are high elevation areas?
#   "The higher-elevation areas tend to be located in the interior of the reserve, especially in the south ... and contain
#   larger amounts of relatively inaccessible forest compared to lower-elevation areas."

    env.data %>% filter(!is.na(elevation_median)) %>%
        ggplot(aes(x = elevation_median, y = distance_to_nature_reserve_boundary)) + geom_point() +
            geom_smooth(method='lm', formula = y ~ x)

    env.data %>% filter(!is.na(elevation_median)) %>%
        mutate(elevation_median.percentile = rank(elevation_median)/n() * 100) %>%
        mutate(elevation_median.quartile = cut(elevation_median.percentile, breaks = c(0,25,50,75,100), labels = c(1,2,3,4))) %>%
        group_by(elevation_median.quartile) %>%
        summarize(distance_to_nature_reserve_boundary.mean = round(mean(distance_to_nature_reserve_boundary) / 10, 0) * 10,
            distance_to_nature_reserve_boundary.sd = round(sd(distance_to_nature_reserve_boundary) / 10, 0) * 10)

    env.data %>% filter(!is.na(elevation_median)) %>%
        lm(distance_to_nature_reserve_boundary ~ elevation_median, data = .) %>% summary()

# save data to file

    save(list="env.data", file=here("rdata","Ailaoshan_environmental.rdata"))
