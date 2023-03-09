library(tidyverse)
mal = read_tsv("allhaps.malathion.200kb.txt.gz")
mal <- na.omit(mal)

myres = mal %>% 
  group_by(chr,pos) %>% 
  nest()
mal2 = mal %>% mutate(treat=str_sub(pool,2,2))
anova1 <- function(df){
  out = anova(lm(asin(sqrt(freq)) ~ treat + founder + treat:founder, data=mal2))
  myF = -pf(out[1,3]/out[2,3],out[1,1],out[2,1],lower.tail=FALSE,
            log.p=TRUE)/log(10)
  myF
}
myres = mal2 %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(logp = map_dbl(data, anova1)) %>%
  select(-data)

anova2 <- function(df){anova(lm(asin(sqrt(freq)) ~ founder + treat %in% founder, data=mal2))}
myres2 = mal %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(logp = map_dbl(data, anova2)) %>%
  select(-data)

anova3 <- function(df){
  out = anova(lm(asin(sqrt(freq)) ~ founder + treat %in% founder, data=mal2))
  myF = -pf(out[1,3]/out[2,3],out[1,1],out[2,1],lower.tail=FALSE,
            log.p=TRUE)/log(10)
  myF
}

myres2 = mal2 %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(logp = map_dbl(data, anova3)) %>%
  select(-data)

join <- full_join(myres, myres2, by = 'pos')
join

