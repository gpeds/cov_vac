#Data treatment

#Packages
pacman::p_load(dplyr,janitor,lubridate,tidyverse,gtsummary, vroom)

df_2022 = vroom("https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SRAG/2022/INFLUD22-03-04-2023.csv") #última atualização opendatasus 	29 de fevereiro de 2024, 12:04 (UTC-03:00)

#Filter dataset
df = clean_names(df_2022)
rm(df_2022)
df = df %>% filter(nu_idade_n < 18)
df1 = df %>% filter(hospital == 1)
df1 = df1 %>% filter(classi_fin==5 & !(cs_sexo == "I"))
df1 = df1 %>% filter(pcr_sars2 == 1|an_sars2==1)
df1 = df1 %>% filter(!(evolucao==3)|!(evolucao==9)) 
df1 = df1 %>% filter(!(suport_ven==9))
df1 = df1 %>% filter(!(uti==9) | !(is.na(uti)))
df1 = df1 %>% filter(!(uti==9))
df1 = df1 %>% filter(!(evolucao==3))
df1 = df1 %>% filter(!(evolucao==9))
df1$asma = df1$asma = NULL
df1$fator_risc = df1$fator_risc = NULL


#Factor to age_class
df1$age <- as.factor(ifelse(df1$nu_idade_n<5,"0-4",
                           ifelse(df1$nu_idade_n<=11.9,"5-11","12-17")))
df1$age <- factor(df1$age, levels = c("0-4", "5-11", "12-17"))

# Recodificando a variável
df1$evolucao <- ifelse(df1$evolucao == 1, "Cura", "Óbito")
df1$uti <- ifelse(df1$uti == 1, "sim", "não")
df1$suport_ven <- ifelse(df1$suport_ven == 1, "TOT", "não")

#Date transformation
df1$dt_sin_pri = dmy(df1$dt_sin_pri)
df1$dt_co_sor = dmy(df1$dt_co_sor)
df1$dt_entuti  = dmy(df1$dt_entuti)
df1$dt_interna  = dmy(df1$dt_interna)
df1$dt_notific  = dmy(df1$dt_notific)
df1$dt_saiduti  = dmy(df1$dt_saiduti)
df1$dt_evoluca  = dmy(df1$dt_evoluca)
df1$dose_1_cov  = dmy(df1$dose_1_cov)
df1$dose_2_cov  = dmy(df1$dose_2_cov)

#Difftime, days to vaccines

df1$vacina_1_dose <- as.numeric(difftime(df1$dt_interna,df1$dose_1_cov, units = "days"))
df1$vacina_2_dose <- as.numeric(difftime(df1$dt_interna,df1$dose_2_cov, units = "days"))
df1$hospdiff <- as.numeric(difftime(df1$dt_interna,df1$dt_sin_pri, units = "days"))
df1$temphosp <- as.numeric(difftime(df1$dt_evoluca,df1$dt_interna, units = "days"))
df1$tempouti <- as.numeric(difftime(df1$dt_saiduti,df1$dt_entuti, units = "days"))



#Variant period criation
df1$variante[df1$dt_notific >= as.Date("2020-03-01") & df1$dt_notific <= as.Date("2021-01-31")] <- "outros"
df1$variante[df1$dt_notific >= as.Date("2021-02-01") & df1$dt_notific <= as.Date("2021-07-31")] <- "gama"
df1$variante[df1$dt_notific >= as.Date("2021-08-01") & df1$dt_notific <= as.Date("2021-12-31")] <- "delta"
df1$variante[df1$dt_notific >= as.Date("2022-01-01") & df1$dt_notific <= as.Date("2023-06-14")] <- "omicron"

#Reorder variant categories
ordem <- c("outros", "gama", "delta", "omicron")
df1$variante <- factor(df1$variante, levels = ordem)
str(df1$variante)


#Numeric variable to variants
df1$variante_numerica <- as.numeric(df1$variante)
summary(factor(df1$variante))

#Comorbidities transformations
comorb_cols <- df1[, 41:53]
df1 <- df1 %>%
  mutate(comor = ifelse(rowSums(comorb_cols == 1, na.rm = TRUE), "Sim",
                              ifelse(rowSums(comorb_cols == 2, na.rm = TRUE), "Não", "Ignorado")))
summary(factor(df1$comor))


#Virus transformations
virus_cols <- df1[, 94:104]

df1 <- df1 %>%
  mutate(coinfec = ifelse(rowSums(virus_cols == 1, na.rm = TRUE), "Sim",
                        ifelse(rowSums(virus_cols < 1, na.rm = TRUE), "Não", "Ignorado")))
summary(factor(df1$coinfec))

#Date transformations
df1$ano = as.integer(format(df1$dt_sin_pri, "%Y"))

#Calcular o número de NA por ano
na_por_ano <- sapply(unique(df1$ano), function(a) sum(is.na(df1$vacina_cov[df1$ano == a])))

#Imprimir o resultado
print(na_por_ano)

#Dataset clean
df1 <- subset(df1, !(is.na(vacina_cov) & ano >= 2021)) #remove linhas do dataframe df1 em com valores NA e o valor da coluna ano >= a 2021.
df1 <- subset(df1, !(vacina_cov == 9 & ano >= 2021)) #remove as linhas do dataframe df1 em com valores 9 e o valor da coluna ano >= a 2021.
na_por_ano <- sapply(unique(df1$ano), function(a) sum(is.na(df1$vacina_cov[df1$ano == a]))) #Calcula o n de valores ausentes NA na coluna vacina_cov de df1 para cada ano único presente na coluna ano e armazena esses valores em um vetor chamado na_por_ano 
df1$vacina_cov[is.na(df1$vacina_cov)] <- 2 #Este comando substitui os valores NA na coluna vacina_cov do dataframe df1 por 2.
summary(factor(df1$vacina_cov))

#Filter
df1 <- df1 %>%
  filter(!(vacina_cov == 1 &  is.na(dose_1_cov)))

df1$infec_vaccine <- as.factor(ifelse(is.na(df1$vacina_1_dose) & is.na(df1$vacina_2_dose),
                                      "não vacinado",
                                      ifelse(df1$vacina_1_dose > 14 & is.na(df1$vacina_2_dose),
                                             "parcialmente vacinado",
                                             ifelse(df1$vacina_1_dose < 14 & is.na(df1$vacina_2_dose),
                                                    "não vacinado",
                                      ifelse(df1$vacina_1_dose > 14 & df1$vacina_2_dose > 14,
                                             "completamente vacinado",
                                             ifelse(df1$vacina_1_dose > 14 & df1$vacina_2_dose < 0,
                                                    "parcialmente vacinado",
                                                    "não vacinado"))))))


df1$infec_vaccine <- relevel(df1$infec_vaccine, ref = "não vacinado")

df1$coluna_concat <- paste0(df1$age, df1$infec_vaccine)

summary(factor(df1$fab_cov_1))

df1 <- df1 %>% 
  mutate(tp_imuno_1 = case_when(
    str_detect(fab_cov_1, "(?i)pfaizer|pizer|fazer|fizer|pfeizer|biotecn|pfizer|pfaizer|pfiser|COMIRNATY|pf|PFIZER COMIRNATY|
               p fizer|biontech|pfazer|pzifer") ~ "pfizer",
    str_detect(fab_cov_1, "(?i)coronavac|corona|butanta|corono|fiocruz|sinovac|butantan|coronavac|coronavac") ~ "coronavac",
    str_detect(fab_cov_1, "(?i)janssen") ~ "janssen",
    is.na(fab_cov_1) ~ "não apresentado",
    TRUE ~ NA_character_
  ))

df1 <- df1 %>% 
  mutate(tp_imuno_2 = case_when(
    str_detect(fab_cov_2, "(?i)pfaizer|pizer|fazer|fizer|pfeizer|biotecn|
               pfizer|pfaizer|pfiser|COMIRNATY|
               p fizer|biontech|pfazer|pzifer") ~ "pfizer",
    str_detect(fab_cov_2, "(?i)coronavac|corona|butanta|corono|fiocruz|sinovac|butantan") ~ "coronavac",
    str_detect(fab_cov_2, "(?i)janssen") ~ "janssen",
    is.na(fab_cov_2) ~ "não apresentado",
    TRUE ~ NA_character_
  ))


df1$coluna_concat <- with(df1, ifelse(infec_vaccine == "não vacinado", "não vacinado",
                                      ifelse(!is.na(tp_imuno_2), paste0(tp_imuno_2, infec_vaccine, age),
                                             ifelse(!is.na(tp_imuno_2), paste0(tp_imuno_1, infec_vaccine),
                                                    NA_character_))))

df1$ano


df1 = df1 %>% filter(!(ano==2020|ano==2021))
df1 = df1 %>%
  filter(!(tp_imuno_1 == "janssen" & !is.na(tp_imuno_1)))
df1 = df1 %>% filter(age=="5-11"|age=="12-17")

df1$regiao <- ifelse(df1$sg_uf %in% c("AM", "RR", "AP", "PA", "TO", "RO", "AC"), "Norte",
                                    ifelse(df1$sg_uf %in% c("MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA"), "Nordeste",
                                           ifelse(df1$sg_uf %in% c("MT", "MS", "GO","DF"), "Centro-Oeste",
                                                  ifelse(df1$sg_uf %in% c("SP", "RJ", "MG", "ES"), "Sudeste",
                                                         ifelse(df1$sg_uf %in% c("PR", "RS", "SC"), "Sul", NA)))))

df1 <- df1 %>%
  mutate(idade_vacina = case_when(
    infec_vaccine == "não vacinado" ~ "não vacinado",
    infec_vaccine == "completamente vacinado" & tp_imuno_2 == "coronavac" ~ "completamente vacinado coronavac",
    infec_vaccine == "completamente vacinado" & tp_imuno_2 == "pfizer" ~ "completamente vacinado pfizer",
    infec_vaccine == "parcialmente vacinado" & tp_imuno_1 == "coronavac" ~ "parcialmente vacinado coronavac",
    infec_vaccine == "parcialmente vacinado" & tp_imuno_1 == "pfizer" ~ "parcialmente vacinado pfizer"))

df1 <- subset(df1, !is.na(idade_vacina))

summary(factor(df1$idade_vacina))
summary(factor(df1$infec_vaccine))


df1$epi_week = epiweek(df1$dt_notific)

attach(df1)

tbl1 = df1 %>% filter(age=="5-11" & between(sem_not,10,34))%>% 
  tbl_summary(include = c(nu_idade_n,
                          cs_sexo,
                          cs_raca,comor,
                          regiao, infec_vaccine, 
                          uti, tempouti, suport_ven,
                          temphosp,
                          tp_imuno_1, tp_imuno_2,
                          coinfec, obesidade,evolucao),
              by=infec_vaccine,
              type = nu_idade_n~"continuous",
              label = list(nu_idade_n~"Idade",
                           cs_sexo~"Sexo",
                           infec_vaccine~"Estrato vacinal",
                           cs_raca~"Raça/Cor",
                           regiao~"Região", 
                           comor~"Comorbidades", 
                           uti~"UTI",
                           tempouti~"Tempo de UTI - dias",
                           suport_ven~"Suporte ventilatório",
                           temphosp~"Tempo de Hospitalização - dias",
                           tp_imuno_1 ~"Tipo de imuno - dose 1",
                           tp_imuno_2~"Tipo de imuno - dose 2",
                           coinfec~"Coinfecção",
                           obesidade~"Obesidade",
                           evolucao~"Óbito"),
              missing = "no") %>% 
  bold_labels() %>% add_p() %>% 
  as_flex_table() %>% 
  flextable::save_as_docx(path="demo_estrato511_simples.docx")

tbl2 = df1 %>% filter(age=="12-17" & between(sem_not,1,22))%>% 
  tbl_summary(include = c(nu_idade_n,
                          cs_sexo,
                          cs_raca,comor,
                          regiao, infec_vaccine, 
                          uti, tempouti, suport_ven,
                          temphosp,
                          coinfec, obesidade,evolucao),
              by=infec_vaccine,
              type = nu_idade_n~"continuous",
              label = list(nu_idade_n~"Idade",
                           cs_sexo~"Sexo",
                           cs_raca~"Raça/Cor",
                           regiao~"Região", 
                           comor~"Comorbidades", 
                           uti~"UTI",
                           tempouti~"Tempo de UTI - dias",
                           suport_ven~"Suporte ventilatório",
                           temphosp~"Tempo de Hospitalização - dias",
                           coinfec~"Coinfecção",
                           obesidade~"Obesidade",
                           evolucao~"Óbito"),
              missing = "no") %>% 
  bold_labels() %>% add_p() %>% 
  as_flex_table() %>% 
  flextable::save_as_docx(path="demo_estrato1217_simples.docx")

tbl3 = df1 %>% filter(age=="5-11"& between(sem_not,10,34))%>% tbl_summary(include = c(nu_idade_n,
                                                               cs_sexo,
                                                               cs_raca,comor,
                                                               regiao, idade_vacina, 
                                                               uti, tempouti, suport_ven,
                                                               temphosp,
                                                               tp_imuno_1, tp_imuno_2,
                                                               coinfec, obesidade,evolucao),
                                                   by=idade_vacina,
                                                   type = nu_idade_n~"continuous",
                                                   label = list(nu_idade_n~"Idade",
                                                                cs_sexo~"Sexo",
                                                                cs_raca~"Raça/Cor",
                                                                regiao~"Região", 
                                                                comor~"Comorbidades", 
                                                                uti~"UTI",
                                                                tempouti~"Tempo de UTI - dias",
                                                                suport_ven~"Suporte ventilatório",
                                                                temphosp~"Tempo de Hospitalização - dias",
                                                                tp_imuno_1 ~"Tipo de imuno - dose 1",
                                                                tp_imuno_2~"Tipo de imuno - dose 2",
                                                                coinfec~"Coinfecção",
                                                                obesidade~"Obesidade",
                                                                evolucao~"Óbito"),
                                                   missing = "no") %>% 
  bold_labels() %>% add_p() %>% 
  as_flex_table() %>% 
  flextable::save_as_docx(path="demo_estrato511_tp.docx")


tbl4 = df1 %>% filter(age=="12-17" & between(sem_not,1,22))%>% tbl_summary(include = c(nu_idade_n,
                                                               cs_sexo,
                                                               cs_raca,comor,
                                                               regiao, idade_vacina, 
                                                               uti, tempouti, suport_ven,
                                                               temphosp,
                                                               coinfec, obesidade,evolucao),
                                                   by=idade_vacina,
                                                   type = nu_idade_n~"continuous",
                                                   label = list(nu_idade_n~"Idade",
                                                                cs_sexo~"Sexo",
                                                                cs_raca~"Raça/Cor",
                                                                regiao~"Região", 
                                                                comor~"Comorbidades", 
                                                                uti~"UTI",
                                                                tempouti~"Tempo de UTI - dias",
                                                                suport_ven~"Suporte ventilatório",
                                                                temphosp~"Tempo de Hospitalização - dias",
                                                                coinfec~"Coinfecção",
                                                                obesidade~"Obesidade",
                                                                evolucao~"Óbito"),
                                                   missing = "no") %>% 
  bold_labels() %>% add_p() %>% 
  as_flex_table() %>% 
  flextable::save_as_docx(path="demo_estrato1217_tp.docx")



#df_epi0 <- df1 %>%
#  group_by(epi_week, idade_vacina) %>%
#  summarise(contagem = n()) %>%
#  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
#  pivot_wider(names_from = idade_vacina, values_from = contagem, values_fill = 0)

#df_epi1 <- df1 %>%
#  group_by(epi_week, infec_vaccine) %>%
#  summarise(contagem = n()) %>%
#  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
#  pivot_wider(names_from = infec_vaccine, values_from = contagem, values_fill = 0)

#df_epi1$`não vacinado` = NULL

#df_epi = cbind(df_epi, df_epi1)
#df_epi$epi_wee.1


#df_0 = readxl::read_xlsx("index.xlsx")
#df_1 = readxl::read_xlsx("pfizer5_11_parcial.xlsx")
#df_2 = readxl::read_xlsx("pfizer5_11_total.xlsx")
#df_3 = readxl::read_xlsx("pfizer12_17_parcial.xlsx")
#df_4 = readxl::read_xlsx("pfizer12_17_total.xlsx")
#df_5 = readxl::read_xlsx("coronavac511_parcial.xlsx")
#df_6 = readxl::read_xlsx("coronavac511_total.xlsx")
#df_7 = readxl::read_xlsx("coronavac1217_parcial.xlsx")
#df_8 = readxl::read_xlsx("coronavac1217_total.xlsx")

#df_pop = left_join(df_0,df_1,by="data")
#df_pop = left_join(df_pop,df_2,by="data")
#df_pop = left_join(df_pop,df_3,by="data")
#df_pop = left_join(df_pop,df_4,by="data")
#df_pop = left_join(df_pop,df_5,by="data")
#df_pop = left_join(df_pop,df_6,by="data")
#df_pop = left_join(df_pop,df_7,by="data")
#df_pop = left_join(df_pop,df_8,by="data")

#df_pop$data = ymd(df_pop$data)
#df_pop$epi_week = epiweek(df_pop$data) 
#df_pop$year = year(df_pop$data) 

#df_pop = df_pop %>% arrange(epi_week,year)

#writexl::write_xlsx(df_pop,"df_pop1.xlsx")

pop = readxl::read_xlsx("df_pop.xlsx")

pop1 <- pop %>%
  select(-c(2:9))

pop = pop %>% arrange(data)


pop$epi_week = epiweek(pop$data) 
pop$year = year(pop$data) 
pop$epi_week[1:2] <- c(1, 1)

pop <- pop %>%
  group_by(epi_week) %>%
  summarise_all(max, na.rm = TRUE)


par(mfrow=c(3,2))
par("mar"=c(4,4, 7,2))

plot(pop$nao_vac_511_cum, main  = "Não vacinados - 5-11 anos", ylab = "Cumulativo", xlab = "", type = "l", lwd=2)
plot(pop$nao_vac_1217_cum, main = "Não vacinados - 12-17 anos", ylab = "", xlab = "", type = "l",lwd=2)
plot(pop$parcial_vac_511_cum, main = "Parcialmente vacinados - 5-11 anos", ylab = "Cumulativo", xlab = "",type = "l", lwd=2)
plot(pop$parcial_vac_1217_cum, main = "Parcialmente vacinados - 12-17 anos", ylab = "", xlab = "",type = "l", lwd=2)
plot(pop$total_vac_511_cum, main = "Totalmente vacinados - 5-11 anos", ylab = "Cumulativo", xlab =  "Semanas epidemiológicas",type = "l", lwd=2)
plot(pop$total_vac_1217_cum, main = "Totalmente vacinados - 12-17 anos", xlab = "Semanas epidemiológicas", ylab="",type = "l", lwd=2)

main_title <- "Soma cumulativa de vacinação por faixa etária e estrato vacinal"
mtext(main_title, side = 3, line = -1.5, outer = TRUE)

par(mfrow=c(2,2))
plot(pop$pfizerparcial_511_cum, main  = "Parcialmente vacinados - 5-11 anos", ylab = "Cumulativo", xlab = "",type = "l", lwd=2)
plot(pop$pfizertotal_511_cum, main = "Totalmente vacinados - 5-11 anos", ylab = "", xlab = "",type = "l", lwd=2)
plot(pop$pfizerparcial_1217_cum, main = "Parcialmente vacinados - 5-11 anos", ylab = "Cumulativo", xlab="Semanas epidemiológicas",type = "l", lwd=2)
plot(pop$pfizertotal_1217_cum, main = "Totalmente vacinados - 12-17 anos", xlab = "Semanas epidemiológicas", ylab="",type = "l", lwd=2)

main_title <- "Soma cumulativa de vacinação com Pfizer e faixa etária"
mtext(main_title, side = 3, line = -1.5, outer = TRUE)

par(mfrow=c(2,2))
plot(pop$coronavacparcial_511_cum, main  = "Parcialmente vacinados - 5-11 anos", ylab = "Cumulativo", xlab = "",type = "l", lwd=2)
plot(pop$coronavactotal_511_cum, main = "Totalmente vacinados - 5-11 anos", ylab = "", xlab = "",type = "l", lwd=2)
plot(pop$coronavacparcial_1217_cum, main = "Parcialmente vacinados - 5-11 anos", ylab = "Cumulativo", xlab="Semanas epidemiológicas",type = "l", lwd=2)
plot(pop$coronavactotal_1217_cum, main = "Totalmente vacinados - 12-17 anos", ylab = "", xlab="Semanas epidemiológicas",type = "l", lwd=2)

main_title <- "Soma cumulativa de vacinação com Coronavac e faixa etária" 
mtext(main_title, side = 3, line = -1.5, outer = TRUE)




#################################
df1$epi_week = epiweek(df1$dt_notific)

df_epi <- df1 %>%
  group_by(epi_week,idade_vacina, age) %>%
  summarise(contagem = n()) %>%
  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
  pivot_wider(names_from = c("idade_vacina","age"), values_from = contagem, values_fill = 0)

df_epi = clean_names(df_epi)

df_epi1 <- df1 %>%
  group_by(epi_week, infec_vaccine, age) %>%
  summarise(contagem = n()) %>%
  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
  pivot_wider(names_from = c("infec_vaccine","age"), values_from = contagem, values_fill = 0)

df_epi1 = clean_names(df_epi1)


df1$epi_week = epiweek(df1$dt_notific)

df_epi2 <- df1 %>%
  group_by(epi_week, idade_vacina, suport_ven, age) %>%
  summarise(contagem = n()) %>%
  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
  pivot_wider(names_from = c("idade_vacina","suport_ven","age"), values_from = contagem, values_fill = 0)

df_epi2$`completamente vacinado pfizer_não_12-17`= NULL
df_epi2$`completamente vacinado pfizer_não_5-11`= NULL
df_epi2$`não vacinado_não_12-17`= NULL
df_epi2$`não vacinado_não_5-11`= NULL
df_epi2$`parcialmente vacinado pfizer_não_5-11`= NULL
df_epi2$`completamente vacinado coronavac_não_12-17`= NULL
df_epi2$`parcialmente vacinado coronavac_não_5-11`= NULL
df_epi2$`parcialmente vacinado coronavac_não_12-17`= NULL
df_epi2$`parcialmente vacinado pfizer_não_12-17`= NULL
df_epi2$`completamente vacinado coronavac_não_5-11`= NULL

df_epi2 = df_epi2 %>% rename("nao_vacinado_tot_511" = `não vacinado_TOT_5-11`,
                             "nao_vacinado_tot_1217" = `não vacinado_TOT_12-17`,
                             "parcial_coronavac_tot_511" = `parcialmente vacinado coronavac_TOT_5-11`,
                             "parcial_coronavac_tot_1217" = `parcialmente vacinado coronavac_TOT_12-17`,
                             "total_corona_tot_511" = `completamente vacinado coronavac_TOT_5-11`,
                             "total_corona_tot_1217" = `completamente vacinado coronavac_TOT_12-17`,
                             "parcial_pfizer_tot_511" = `parcialmente vacinado pfizer_TOT_5-11`,
                             "parcial_pfizer_tot_1217" = `parcialmente vacinado pfizer_TOT_12-17`,
                             "total_pfizer_tot_511" = `completamente vacinado pfizer_TOT_5-11`,
                             "total_pfizer_tot_1217" = `completamente vacinado pfizer_TOT_12-17`,
                             )

df_epi3 <- df1 %>%
  group_by(epi_week, infec_vaccine, suport_ven, age) %>%
  summarise(contagem = n()) %>%
  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
  pivot_wider(names_from = c("infec_vaccine","suport_ven","age"), values_from = contagem, values_fill = 0)

df_epi3$`não vacinado_não_12-17`= NULL
df_epi3$`não vacinado_não_5-11`= NULL
df_epi3$`completamente vacinado_não_5-11`= NULL
df_epi3$`completamente vacinado_não_12-17`= NULL
df_epi3$`parcialmente vacinado_não_12-17`= NULL
df_epi3$`parcialmente vacinado_não_5-11`= NULL  

  
df_epi3 = df_epi3 %>% rename("nao_vacinado_tot_511" = `não vacinado_TOT_5-11`,
                             "nao_vacinado_tot_1217" = `não vacinado_TOT_12-17`,
                             "parcial_tot_511" = `parcialmente vacinado_TOT_5-11`,
                             "parcial_tot_1217" = `parcialmente vacinado_TOT_12-17`,
                             "total_tot_1217" = `completamente vacinado_TOT_12-17`,
                             "total_tot_511" = `completamente vacinado_TOT_5-11`)



df_epi4 <- df1 %>%
  group_by(epi_week, idade_vacina, evolucao, age) %>%
  summarise(contagem = n()) %>%
  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
  pivot_wider(names_from = c("idade_vacina","evolucao","age"), values_from = contagem, values_fill = 0)

df_epi4 = clean_names(df_epi4)

df_epi4$completamente_vacinado_pfizer_cura_12_17 = NULL
df_epi4$nao_vacinado_cura_5_11 = NULL
df_epi4$nao_vacinado_cura_12_17= NULL
df_epi4$parcialmente_vacinado_coronavac_cura_5_11= NULL
df_epi4$parcialmente_vacinado_pfizer_cura_12_17= NULL
df_epi4$parcialmente_vacinado_pfizer_cura_5_11= NULL
df_epi4$parcialmente_vacinado_coronavac_cura_12_17= NULL
df_epi4$parcialmente_vacinado_pfizer_cura_12_17= NULL
df_epi4$completamente_vacinado_coronavac_cura_12_17= NULL
df_epi4$completamente_vacinado_coronavac_cura_5_11= NULL


df_epi5 <- df1 %>%
  group_by(epi_week, infec_vaccine, evolucao, age) %>%
  summarise(contagem = n()) %>%
  ungroup() %>%
  # Pivotar os dados para que cada tipo de vacina seja uma coluna
  pivot_wider(names_from = c("infec_vaccine","evolucao","age"), values_from = contagem, values_fill = 0)

df_epi5 = clean_names(df_epi5)

df_epi5$nao_vacinado_cura_5_11 = NULL
df_epi5$nao_vacinado_cura_12_17 = NULL
df_epi5$completamente_vacinado_cura_12_17 = NULL
df_epi5$completamente_vacinado_cura_5_11 = NULL
df_epi5$parcialmente_vacinado_cura_12_17 = NULL
df_epi5$parcialmente_vacinado_cura_5_11 = NULL

#------------------------------------------------------------------------------------------------------------------------------------#

df_stata = left_join(pop, df_epi, by="epi_week")
df_stata = left_join(df_stata, df_epi2, by="epi_week")
df_stata = left_join(df_stata, df_epi4, by="epi_week")

df_stata$taxa_pfpr_511_srag = df_stata$parcialmente_vacinado_pfizer_5_11/df_stata$pfizerparcial_511_cum*100000
df_stata$taxa_pfpr_1217_srag = df_stata$parcialmente_vacinado_pfizer_12_17/df_stata$pfizerparcial_1217_cum*100000
df_stata$taxa_pftt_511_srag = df_stata$completamente_vacinado_pfizer_5_11/df_stata$pfizertotal_511_cum*100000
df_stata$taxa_pftt_1217_srag = df_stata$completamente_vacinado_pfizer_12_17/df_stata$pfizertotal_1217_cum*100000
df_stata$taxa_cvpr_511_srag = df_stata$parcialmente_vacinado_coronavac_5_11/df_stata$pfizerparcial_511_cum*100000
df_stata$taxa_cvpr_1217_srag = df_stata$parcialmente_vacinado_coronavac_12_17/df_stata$pfizerparcial_1217_cum*100000
df_stata$taxa_cvtt_511_srag = df_stata$completamente_vacinado_coronavac_5_11/df_stata$pfizertotal_511_cum*100000
df_stata$taxa_cvtt_1217_srag = df_stata$completamente_vacinado_coronavac_12_17/df_stata$pfizertotal_1217_cum*100000
df_stata$taxa_nv_511_srag = df_stata$nao_vacinado_5_11/df_stata$nao_vac_511_cum*100000
df_stata$taxa_nv_1217_srag = df_stata$nao_vacinado_12_17/df_stata$nao_vac_1217_cum*100000

df_stata$taxa_pfpr_511_tot = df_stata$parcial_pfizer_tot_511/df_stata$pfizerparcial_511_cum*100000
df_stata$taxa_pfpr_1217_tot = df_stata$parcial_pfizer_tot_1217/df_stata$pfizerparcial_1217_cum*100000
df_stata$taxa_pftt_511_tot = df_stata$total_pfizer_tot_511/df_stata$pfizertotal_511_cum*100000
df_stata$taxa_pftt_1217_tot = df_stata$total_pfizer_tot_1217/df_stata$pfizertotal_1217_cum*100000
df_stata$taxa_cvpr_511_tot = df_stata$parcial_coronavac_tot_511/df_stata$coronavacparcial_511_cum*100000
df_stata$taxa_cvpr_1217_tot = df_stata$parcial_coronavac_tot_1217/df_stata$coronavacparcial_1217_cum*100000
df_stata$taxa_cvtt_511_tot = df_stata$total_corona_tot_511/df_stata$coronavactotal_511_cum*100000
df_stata$taxa_cvtt_1217_tot = df_stata$total_corona_tot_1217/df_stata$coronavactotal_1217_cum*100000
df_stata$taxa_nv_511_tot = df_stata$nao_vacinado_tot_511/df_stata$nao_vac_511_cum*100000
df_stata$taxa_nv_1217_tot = df_stata$nao_vacinado_tot_1217/df_stata$nao_vac_1217_cum*100000

df_stata$taxa_pfpr_511_obito = df_stata$parcialmente_vacinado_pfizer_obito_5_11/df_stata$pfizerparcial_511_cum*100000
df_stata$taxa_pfpr_1217_obito = df_stata$parcialmente_vacinado_pfizer_obito_12_17/df_stata$pfizerparcial_1217_cum*100000
df_stata$taxa_pftt_511_obito = df_stata$completamente_vacinado_pfizer_obito_5_11/df_stata$pfizertotal_511_cum*100000
df_stata$taxa_pftt_1217_obito = df_stata$completamente_vacinado_pfizer_obito_12_17/df_stata$pfizertotal_1217_cum*100000
df_stata$taxa_cvpr_511_obito = df_stata$parcialmente_vacinado_coronavac_obito_5_11/df_stata$coronavacparcial_511_cum*100000
df_stata$taxa_cvpr_1217_obito = df_stata$parcialmente_vacinado_coronavac_obito_12_17/df_stata$coronavacparcial_1217_cum*100000
df_stata$taxa_cvtt_511_obito = df_stata$completamente_vacinado_coronavac_obito_5_11/df_stata$coronavactotal_511_cum*100000
df_stata$taxa_cvtt_1217_obito = df_stata$completamente_vacinado_coronavac_obito_12_17/df_stata$coronavactotal_1217_cum*100000
df_stata$taxa_nv_511_obito = df_stata$nao_vacinado_obito_5_11/df_stata$nao_vac_511_cum*100000
df_stata$taxa_nv_1217_obito = df_stata$nao_vacinado_obito_12_17/df_stata$nao_vac_1217_cum*100000



###############
df_stata1 = left_join(pop, df_epi1, by="epi_week")
df_stata1 = left_join(df_stata1, df_epi3, by="epi_week")
df_stata1 = left_join(df_stata1, df_epi5, by="epi_week")

df_stata1$taxa_pr_511_srag = df_stata1$parcialmente_vacinado_5_11/df_stata1$parcial_vac_511_cum*100000
df_stata1$taxa_pr_1217_srag = df_stata1$parcialmente_vacinado_12_17/df_stata1$parcial_vac_1217_cum*100000
df_stata1$taxa_tt_511_srag = df_stata1$completamente_vacinado_5_11/df_stata1$total_vac_511_cum*100000
df_stata1$taxa_tt_1217_srag = df_stata1$completamente_vacinado_12_17/df_stata1$total_vac_1217_cum*100000

df_stata1$taxa_pr_511_tot = df_stata1$parcial_tot_511/df_stata1$parcial_vac_511_cum*100000
df_stata1$taxa_pr_1217_tot = df_stata1$parcial_tot_1217/df_stata1$parcial_vac_1217_cum*100000
df_stata1$taxa_tt_511_tot = df_stata1$total_tot_511/df_stata1$total_vac_511_cum*100000
df_stata1$taxa_tt_1217_tot = df_stata1$total_tot_1217/df_stata1$total_vac_1217_cum*100000


df_stata1$taxa_pr_511_obito = df_stata1$parcialmente_vacinado_obito_5_11/df_stata1$parcial_vac_511_cum*100000
df_stata1$taxa_pr_1217_obito = df_stata1$parcialmente_vacinado_obito_12_17/df_stata1$parcial_vac_1217_cum*100000
df_stata1$taxa_tt_511_obito = df_stata1$completamente_vacinado_obito_5_11/df_stata1$total_vac_511_cum*100000
df_stata1$taxa_tt_1217_obito = df_stata1$completamente_vacinado_obito_12_17/df_stata1$total_vac_1217_cum*100000

df_stata1$taxa_nv_511_srag = df_stata1$nao_vacinado_5_11/df_stata1$nao_vac_511_cum*100000
df_stata1$taxa_nv_1217_srag = df_stata1$nao_vacinado_12_17/df_stata1$nao_vac_1217_cum*100000
df_stata1$taxa_nv_511_tot = df_stata1$nao_vacinado_tot_511/df_stata1$nao_vac_511_cum*100000
df_stata1$taxa_nv_1217_tot = df_stata1$nao_vacinado_tot_1217/df_stata1$nao_vac_1217_cum*100000
df_stata1$taxa_nv_511_obito = df_stata1$nao_vacinado_obito_5_11/df_stata1$nao_vac_511_cum*100000
df_stata1$taxa_nv_1217_obito = df_stata1$nao_vacinado_obito_12_17/df_stata1$nao_vac_1217_cum*100000


##############
df_511 = df_stata[10:34,]
df_1217 = df_stata[1:22,]

options(scipen = 999)

srag_511 <- df_511 %>% select(epi_week, taxa_pfpr_511_srag, taxa_pftt_511_srag, taxa_nv_511_srag, taxa_cvpr_511_srag, taxa_cvtt_511_srag) %>% 
  pivot_longer(
    cols = c(taxa_pfpr_511_srag, taxa_pftt_511_srag, taxa_nv_511_srag, taxa_cvpr_511_srag, taxa_cvtt_511_srag),
    names_to = "Vacinacao",
    values_to = "valor"
  )

srag_511$Vacinacao = factor(srag_511$Vacinacao)
srag_511$Vacinacao = relevel(srag_511$Vacinacao, ref = "taxa_nv_511_srag")

model_511_srag = glm(valor~Vacinacao, family = "quasipoisson", data=srag_511)


# Calcular os odds ratio
irr_511_srag <- exp(coef(model_511_srag))

confint_511_srag <- exp(confint(model_511_srag))


# Obter os valores de p
p_value_511_srag <- c("0.48","0.01","0.85","0.01","0.30")

# Criar um dataframe com os resultados
resultados_511_srag <- data.frame(
  relative_risk = irr_511_srag,
  p_value = p_value_511_srag,  # Adicionando o valor de p
  lower_ci = confint_511_srag[, 1],
  upper_ci = confint_511_srag[, 2]
)

print(resultados_511_srag)



###############
tot_511 <- df_511 %>% select(epi_week, taxa_pfpr_511_tot, taxa_pftt_511_tot, taxa_nv_511_tot, taxa_cvpr_511_tot, taxa_cvtt_511_tot) %>% 
  pivot_longer(
    cols = c(taxa_pfpr_511_tot, taxa_pftt_511_tot, taxa_nv_511_tot, taxa_cvpr_511_tot, taxa_cvtt_511_tot),
    names_to = "Vacinacao",
    values_to = "valor"
  )

tot_511$Vacinacao = factor(tot_511$Vacinacao)
tot_511$Vacinacao = relevel(tot_511$Vacinacao, ref = "taxa_nv_511_tot")

model_511_tot = glm(valor~Vacinacao, family = "quasipoisson", data=tot_511)


# Calcular os odds ratio
irr_511_tot <- exp(coef(model_511_tot))

confint_511_tot <- exp(confint(model_511_tot))


# Obter os valores de p
p_value_511_tot <- c("<0.01","<0.01","<0.01","<0.01","<0.01")

# Criar um dataframe com os resultados
resultados_511_tot <- data.frame(
  relative_risk = irr_511_tot,
  p_value = p_value_511_tot,  # Adicionando o valor de p
  lower_ci = confint_511_tot[, 1],
  upper_ci = confint_511_tot[, 2]
)

print(resultados_511_tot)

########################

obito_511 <- df_511 %>% select(epi_week, taxa_pfpr_511_obito, taxa_pftt_511_obito, taxa_nv_511_obito, taxa_cvpr_511_obito, taxa_cvtt_511_obito) %>% 
  pivot_longer(
    cols = c(taxa_pfpr_511_obito, taxa_pftt_511_obito, taxa_nv_511_obito, taxa_cvpr_511_obito, taxa_cvtt_511_obito),
    names_to = "Vacinacao",
    values_to = "valor"
  )

obito_511$Vacinacao = factor(obito_511$Vacinacao)
obito_511$Vacinacao = relevel(obito_511$Vacinacao, ref = "taxa_nv_511_obito")

model_511_obito = glm(valor~Vacinacao, family = "quasipoisson", data=obito_511)


# Calcular os odds ratio
irr_511_obito <- exp(coef(model_511_obito))

confint_511_obito <- exp(confint(model_511_obito))


# Obter os valores de p
p_value_511_obito <- c("<0.01","0.99","0.99","<0.01","<0.01")

# Criar um dataframe com os resultados
resultados_511_obito <- data.frame(
  relative_risk = irr_511_obito,
  p_value = p_value_511_obito,  # Adicionando o valor de p
  lower_ci = NA,
  upper_ci = NA
)

print(resultados_511_obito)

dados_511 = rbind(resultados_511_srag,resultados_511_tot,resultados_511_obito)

##################################################

srag_1217 <- df_1217 %>% select(epi_week, taxa_pfpr_1217_srag, taxa_pftt_1217_srag, taxa_nv_1217_srag, taxa_cvpr_1217_srag, taxa_cvtt_1217_srag) %>% 
  pivot_longer(
    cols = c(taxa_pfpr_1217_srag, taxa_pftt_1217_srag, taxa_nv_1217_srag, taxa_cvpr_1217_srag, taxa_cvtt_1217_srag),
    names_to = "Vacinacao",
    values_to = "valor"
  )

srag_1217$Vacinacao = factor(srag_1217$Vacinacao)
srag_1217$Vacinacao = relevel(srag_1217$Vacinacao, ref = "taxa_nv_1217_srag")

model_1217_srag = glm(valor~Vacinacao, family = "quasipoisson", data=srag_1217)


# Calcular os odds ratio
irr_1217_srag <- exp(coef(model_1217_srag))

confint_1217_srag <- exp(confint(model_1217_srag))


# Obter os valores de p
p_value_1217_srag <- c("<0.01","<0.01","<0.01","<0.07","<0.01")

# Criar um dataframe com os resultados
resultados_1217_srag <- data.frame(
  relative_risk = irr_1217_srag,
  p_value = p_value_1217_srag,  # Adicionando o valor de p
  lower_ci = confint_1217_srag[, 1],
  upper_ci = confint_1217_srag[, 2]
)

print(resultados_1217_srag)



###############
tot_1217 <- df_1217 %>% select(epi_week, taxa_pfpr_1217_tot, taxa_pftt_1217_tot, taxa_nv_1217_tot, taxa_cvpr_1217_tot, taxa_cvtt_1217_tot) %>% 
  pivot_longer(
    cols = c(taxa_pfpr_1217_tot, taxa_pftt_1217_tot, taxa_nv_1217_tot, taxa_cvpr_1217_tot, taxa_cvtt_1217_tot),
    names_to = "Vacinacao",
    values_to = "valor"
  )

tot_1217$Vacinacao = factor(tot_1217$Vacinacao)
tot_1217$Vacinacao = relevel(tot_1217$Vacinacao, ref = "taxa_nv_1217_tot")

model_1217_tot = glm(valor~Vacinacao, family = "quasipoisson", data=tot_1217)


# Calcular os odds ratio
irr_1217_tot <- exp(coef(model_1217_tot))

confint_1217_tot <- exp(confint(model_1217_tot))


# Obter os valores de p
p_value_1217_tot <- c("<0.01","0.99","0.99","0.05","<0.01")

# Criar um dataframe com os resultados
resultados_1217_tot <- data.frame(
  relative_risk = irr_1217_tot,
  p_value = p_value_1217_tot,  # Adicionando o valor de p
  lower_ci = NA,
  upper_ci = NA
)

print(resultados_1217_tot)

########################

obito_1217 <- df_1217 %>% select(epi_week, taxa_pfpr_1217_obito, taxa_pftt_1217_obito, taxa_nv_1217_obito, taxa_cvpr_1217_obito, taxa_cvtt_1217_obito) %>% 
  pivot_longer(
    cols = c(taxa_pfpr_1217_obito, taxa_pftt_1217_obito, taxa_nv_1217_obito, taxa_cvpr_1217_obito, taxa_cvtt_1217_obito),
    names_to = "Vacinacao",
    values_to = "valor"
  )

obito_1217$Vacinacao = factor(obito_1217$Vacinacao)
obito_1217$Vacinacao = relevel(obito_1217$Vacinacao, ref = "taxa_nv_1217_obito")

model_1217_obito = glm(valor~Vacinacao, family = "quasipoisson", data=obito_1217)


# Calcular os odds ratio
irr_1217_obito <- exp(coef(model_1217_obito))

confint_1217_obito <- exp(confint(model_1217_obito))


# Obter os valores de p
p_value_1217_obito <- c("<0.01","0.99","0.99","<0.20","<0.17")

# Criar um dataframe com os resultados
resultados_1217_obito <- data.frame(
  relative_risk = irr_1217_obito,
  p_value = p_value_1217_obito,  # Adicionando o valor de p
  lower_ci = NA,
  upper_ci = NA
)

print(resultados_1217_obito)

dados_1217 = rbind(resultados_1217_srag,resultados_1217_tot,resultados_1217_obito)


########################


df_511_b = df_stata1[10:34,]
df_1217_b = df_stata1[1:22,]

options(scipen = 999)

srag_511b <- df_511_b %>% select(epi_week, taxa_pr_511_srag, taxa_tt_511_srag, taxa_nv_511_srag) %>% 
  pivot_longer(
    cols = c(taxa_pr_511_srag, taxa_tt_511_srag, taxa_nv_511_srag),
    names_to = "Vacinacao",
    values_to = "valor"
  )

srag_511b$Vacinacao = factor(srag_511b$Vacinacao)
srag_511b$Vacinacao = relevel(srag_511b$Vacinacao, ref = "taxa_nv_511_srag")

model_511_sragb = glm(valor~Vacinacao, family = "quasipoisson", data=srag_511b)

# Calcular os odds ratio
irr_511_sragb <- exp(coef(model_511_sragb))

confint_511_sragb <- exp(confint(model_511_sragb))


# Obter os valores de p
p_value_511_sragb <- c("0.02","<0.01","<0.01")

# Criar um dataframe com os resultados
resultados_511_sragb <- data.frame(
  relative_risk = irr_511_sragb,
  p_value = p_value_511_sragb,  # Adicionando o valor de p
  lower_ci = confint_511_sragb[, 1],
  upper_ci = confint_511_sragb[, 2]
)

print(resultados_511_sragb)



###############
tot_511b <- df_511_b %>% select(epi_week, taxa_pr_511_tot, taxa_tt_511_tot, taxa_nv_511_tot) %>% 
  pivot_longer(
    cols = c(taxa_pr_511_tot, taxa_tt_511_tot, taxa_nv_511_tot),
    names_to = "Vacinacao",
    values_to = "valor"
  )

tot_511b$Vacinacao = factor(tot_511b$Vacinacao)
tot_511b$Vacinacao = relevel(tot_511b$Vacinacao, ref = "taxa_nv_511_tot")

model_511_totb = glm(valor~Vacinacao, family = "quasipoisson", data=tot_511b)


# Calcular os odds ratio
irr_511_totb <- exp(coef(model_511_totb))

confint_511_totb <- exp(confint(model_511_totb))


# Obter os valores de p
p_value_511_totb <- c("<0.01","<0.01","<0.01")

# Criar um dataframe com os resultados
resultados_511_totb <- data.frame(
  relative_risk = irr_511_totb,
  p_value = p_value_511_totb,  # Adicionando o valor de p
  lower_ci = confint_511_totb[, 1],
  upper_ci = confint_511_totb[, 2])

print(resultados_511_totb)

########################

obito_511b <- df_511_b %>% select(epi_week, taxa_pr_511_obito, taxa_tt_511_obito, taxa_nv_511_obito) %>% 
  pivot_longer(
    cols = c(taxa_pr_511_obito, taxa_tt_511_obito, taxa_nv_511_obito),
    names_to = "Vacinacao",
    values_to = "valor"
  )

obito_511b$Vacinacao = factor(obito_511b$Vacinacao)
obito_511b$Vacinacao = relevel(obito_511b$Vacinacao, ref = "taxa_nv_511_obito")

model_511_obitob = glm(valor~Vacinacao, family = "quasipoisson", data=obito_511b)


# Calcular os odds ratio
irr_511_obitob <- exp(coef(model_511_obitob))

confint_511_obitob <- exp(confint(model_511_obitob))


# Obter os valores de p
p_value_511_obitob <- c("<0.01","<0.01","<0.01")

# Criar um dataframe com os resultados
resultados_511_obitob <- data.frame(
  relative_risk = irr_511_obitob,
  p_value = p_value_511_obitob,  # Adicionando o valor de p
  lower_ci = confint_511_obitob[,1],
  upper_ci = confint_511_obitob[,2]
)

print(resultados_511_obitob)

dados_511b = rbind(resultados_511_sragb,resultados_511_totb,resultados_511_obitob)

##################################################

srag_1217b <- df_1217_b %>% select(epi_week, taxa_pr_1217_srag, taxa_tt_1217_srag, taxa_nv_1217_srag) %>% 
  pivot_longer(
    cols = c(taxa_pr_1217_srag, taxa_tt_1217_srag, taxa_nv_1217_srag),
    names_to = "Vacinacao",
    values_to = "valor"
  )

srag_1217b$Vacinacao = factor(srag_1217b$Vacinacao)
srag_1217b$Vacinacao = relevel(srag_1217b$Vacinacao, ref = "taxa_nv_1217_srag")

model_1217_sragb = glm(valor~Vacinacao, family = "quasipoisson", data=srag_1217b)

# Calcular os odds ratio
irr_1217_sragb <- exp(coef(model_1217_sragb))

confint_1217_sragb <- exp(confint(model_1217_sragb))


# Obter os valores de p
p_value_1217_sragb <- c("<0.01","<0.01","<0.01")

# Criar um dataframe com os resultados
resultados_1217_sragb <- data.frame(
  relative_risk = irr_1217_sragb,
  p_value = p_value_1217_sragb,  # Adicionando o valor de p
  lower_ci = confint_1217_sragb[, 1],
  upper_ci = confint_1217_sragb[, 2]
)

print(resultados_1217_sragb)



###############
tot_1217b <- df_1217_b %>% select(epi_week, taxa_pr_1217_tot, taxa_tt_1217_tot, taxa_nv_1217_tot) %>% 
  pivot_longer(
    cols = c(taxa_pr_1217_tot, taxa_tt_1217_tot, taxa_nv_1217_tot),
    names_to = "Vacinacao",
    values_to = "valor"
  )

tot_1217b$Vacinacao = factor(tot_1217b$Vacinacao)
tot_1217b$Vacinacao = relevel(tot_1217b$Vacinacao, ref = "taxa_nv_1217_tot")

model_1217_totb = glm(valor~Vacinacao, family = "quasipoisson", data=tot_1217b)


# Calcular os odds ratio
irr_1217_totb <- exp(coef(model_1217_totb))

confint_1217_totb <- exp(confint(model_1217_totb))


# Obter os valores de p
p_value_1217_totb <- c("<0.01","0.02","0.01")

# Criar um dataframe com os resultados
resultados_1217_totb <- data.frame(
  relative_risk = irr_1217_totb,
  p_value = p_value_1217_totb,  # Adicionando o valor de p
  lower_ci = confint_1217_totb[, 1],
  upper_ci = confint_1217_totb[, 2])

print(resultados_1217_totb)

########################

obito_1217b <- df_1217_b %>% select(epi_week, taxa_pr_1217_obito, taxa_tt_1217_obito, taxa_nv_1217_obito) %>% 
  pivot_longer(
    cols = c(taxa_pr_1217_obito, taxa_tt_1217_obito, taxa_nv_1217_obito),
    names_to = "Vacinacao",
    values_to = "valor"
  )

obito_1217b$Vacinacao = factor(obito_1217b$Vacinacao)
obito_1217b$Vacinacao = relevel(obito_1217b$Vacinacao, ref = "taxa_nv_1217_obito")

model_1217_obitob = glm(valor~Vacinacao, family = "quasipoisson", data=obito_1217b)


# Calcular os odds ratio
irr_1217_obitob <- exp(coef(model_1217_obitob))

confint_1217_obitob <- exp(confint(model_1217_obitob))


# Obter os valores de p
p_value_1217_obitob <- c("<0.01","<0.01","<0.01")

# Criar um dataframe com os resultados
resultados_1217_obitob <- data.frame(
  relative_risk = irr_1217_obitob,
  p_value = p_value_1217_obitob,  # Adicionando o valor de p
  lower_ci = confint_1217_obitob[,1],
  upper_ci = confint_1217_obitob[,2]
)

print(resultados_1217_obitob)

dados_1217b = rbind(resultados_1217_sragb,resultados_1217_totb,resultados_1217_obitob)

openxlsx::write.xlsx(dados_511, "dados_511_tpimuno.xlsx", rowNames=TRUE)
openxlsx::write.xlsx(dados_511b, "dados_511_geral.xlsx", rowNames=TRUE)
openxlsx::write.xlsx(dados_1217, "dados_1217_tpimuno.xlsx", rowNames=TRUE)
openxlsx::write.xlsx(dados_1217b, "dados_1217_geral.xlsx", rowNames=TRUE)

########################
srag_511$Vacinacao <- factor(srag_511$Vacinacao, levels = c("taxa_nv_511_srag", "taxa_cvpr_511_srag", "taxa_cvtt_511_srag", 
                                                      "taxa_pfpr_511_srag", "taxa_pftt_511_srag"), 
                          labels = c("Não vacinados", "Parcialmente - coronavac", "Totalmente - coronavac", "Parcialmente - pfizer", 
                                     "Totalmente - pfizer"))


ggplot(srag_511, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de SRAG por tipo dose e tipo de imunizante - 5-11 anos") +
  theme_minimal()

srag_1217$Vacinacao <- factor(srag_1217$Vacinacao, levels = c("taxa_nv_1217_srag", "taxa_cvpr_1217_srag", "taxa_cvtt_1217_srag", 
                                                            "taxa_pfpr_1217_srag", "taxa_pftt_1217_srag"), 
                             labels = c("Não vacinados", "Parcialmente - coronavac", "Totalmente - coronavac", "Parcialmente - pfizer", 
                                        "Totalmente - pfizer"))


ggplot(srag_1217, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de SRAG por tipo dose e tipo de imunizante - 12-17 anos") +
  theme_minimal()

########################################################################################


tot_511$Vacinacao <- factor(tot_511$Vacinacao, levels = c("taxa_nv_511_tot", "taxa_cvpr_511_tot", "taxa_cvtt_511_tot", 
                                                            "taxa_pfpr_511_tot", "taxa_pftt_511_tot"), 
                             labels = c("Não vacinados", "Parcialmente - coronavac", "Totalmente - coronavac", "Parcialmente - pfizer", 
                                        "Totalmente - pfizer"))


ggplot(tot_511, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de tot por tipo dose e tipo de imunizante - 5-11 anos") +
  theme_minimal()

tot_1217$Vacinacao <- factor(tot_1217$Vacinacao, levels = c("taxa_nv_1217_tot", "taxa_cvpr_1217_tot", "taxa_cvtt_1217_tot", 
                                                              "taxa_pfpr_1217_tot", "taxa_pftt_1217_tot"), 
                              labels = c("Não vacinados", "Parcialmente - coronavac", "Totalmente - coronavac", "Parcialmente - pfizer", 
                                         "Totalmente - pfizer"))


ggplot(tot_1217, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de tot por tipo dose e tipo de imunizante - 12-17 anos") +
  theme_minimal()

#####################################

obito_511$Vacinacao <- factor(obito_511$Vacinacao, levels = c("taxa_nv_511_obito", "taxa_cvpr_511_obito", "taxa_cvtt_511_obito", 
                                                          "taxa_pfpr_511_obito", "taxa_pftt_511_obito"), 
                            labels = c("Não vacinados", "Parcialmente - coronavac", "Totalmente - coronavac", "Parcialmente - pfizer", 
                                       "Totalmente - pfizer"))


ggplot(obito_511, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de tot por tipo dose e tipo de imunizante - 5-11 anos") +
  theme_minimal()

obito_1217$Vacinacao <- factor(obito_1217$Vacinacao, levels = c("taxa_nv_1217_obito", "taxa_cvpr_1217_obito", "taxa_cvtt_1217_obito", 
                                                            "taxa_pfpr_1217_obito", "taxa_pftt_1217_obito"), 
                             labels = c("Não vacinados", "Parcialmente - coronavac", "Totalmente - coronavac", "Parcialmente - pfizer", 
                                        "Totalmente - pfizer"))


ggplot(obito_1217, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de óbito por tipo dose e tipo de imunizante - 12-17 anos") +
  theme_minimal()



########################
srag_511b$Vacinacao <- factor(srag_511b$Vacinacao, levels = c("taxa_pr_511_srag", "taxa_tt_511_srag", "taxa_nv_511_srag"), 
                             labels = c("Parcialmente vacinado", "Totalmente vacinado", "Não vacinado"))


ggplot(srag_511b, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de SRAG por estrato vacinal - 5-11 anos") +
  theme_minimal()

srag_1217b$Vacinacao <- factor(srag_1217b$Vacinacao, levels = c("taxa_pr_1217_srag", "taxa_tt_1217_srag", "taxa_nv_1217_srag"), 
                              labels = c("Parcialmente vacinado", "Totalmente vacinado", "Não vacinado"))


ggplot(srag_1217b, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de SRAG por estrato vacinal - 12-17 anos") +
  theme_minimal()


######


tot_511b$Vacinacao <- factor(tot_511b$Vacinacao, levels = c("taxa_pr_511_tot", "taxa_tt_511_tot", "taxa_nv_511_tot"), 
                              labels = c("Parcialmente vacinado", "Totalmente vacinado", "Não vacinado"))


ggplot(tot_511b, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de tot por estrato vacinal - 5-11 anos") +
  theme_minimal()

tot_1217b$Vacinacao <- factor(tot_1217b$Vacinacao, levels = c("taxa_pr_1217_tot", "taxa_tt_1217_tot", "taxa_nv_1217_tot"), 
                               labels = c("Parcialmente vacinado", "Totalmente vacinado", "Não vacinado"))


ggplot(tot_1217b, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de tot por estrato vacinal - 12-17 anos") +
  theme_minimal()

####

obito_511b$Vacinacao <- factor(obito_511b$Vacinacao, levels = c("taxa_pr_511_obito", "taxa_tt_511_obito", "taxa_nv_511_obito"), 
                              labels = c("Parcialmente vacinado", "Totalmente vacinado", "Não vacinado"))


ggplot(obito_511b, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de óbito por estrato vacinal - 5-11 anos") +
  theme_minimal()

obito_1217b$Vacinacao <- factor(obito_1217b$Vacinacao, levels = c("taxa_pr_1217_obito", "taxa_tt_1217_obito", "taxa_nv_1217_obito"), 
                               labels = c("Parcialmente vacinado", "Totalmente vacinado", "Não vacinado"))


ggplot(obito_1217b, aes(x = epi_week, y = valor, color = Vacinacao)) +
  geom_line() +
  labs(x = "Semana Epidemiológica", y = "Valor", color = "Tipo de Vacinação",
       title = "Taxa de óbito por estrato vacinal - 12-17 anos") +
  theme_minimal()
