#######################
## Recapture rates ####
#######################


test<-data.frame()
a<- 1

for(a in 1: length(data_list) ){

  fragments <- unique(capture_data$site)
  
  
id<-data_list[[a]]$ids    
captures<-  sum(data_list[[a]]$CH)


n.sesions<- data_list[[a]]$n.sessions
fragment<-data_list[[a]]$fragment
sex<-data_list[[a]]$sex

id.data<-data.frame(sex, fragment, n.sesions, captures, id)
test<-rbind(test, id.data)

}


test %>% summarise(mean(captures))

test %>% group_by(sex) %>% summarise(mean(captures))

table.cap<-test %>%   group_by(fragment, sex) %>% 
  summarise(mean(captures)) 
  
table.cap <- table.cap %>%  mutate(sex= ifelse(sex== "1", "female", "male" )) %>% 
                            mutate(protected= ifelse(fragment == "M13" | fragment == "M20", "disturbed", "protected" )) 

write.csv(table.cap, file =  file.path(analysisDir, modelName, "Recaptures_means_grid.csv"))


print(xtable(table.cap, type = "latex",
             align = paste(c("l",rep("c",ncol(table.cap))),collapse = "")),
      floating = FALSE, include.colnames = FALSE, include.rownames = TRUE,
      add.to.row = addtorow,
      file = file.path(analysisDir, modelName, "Recaptures_means_grid.tex"))
