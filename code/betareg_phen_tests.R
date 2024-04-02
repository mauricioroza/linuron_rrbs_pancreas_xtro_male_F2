library(betareg)

perc <- df %>% mutate_if(is.numeric, list(~ . / 100))

t18.1 <- betareg(`18:1` ~ pnliprp2 + lmf1, data = perc)

summary(t18.1)

library(ggplot2)
ggplot(perc, aes(x = vti1a, y = `18:1`)) +
  geom_point(size = 4, aes(fill = treatment), shape = 21) +
  scale_fill_grey() +
  geom_line(aes(y = predict(t18.1, perc), 
                colour = "logit", linetype = "logit")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()

t16.1 <- betareg(`16:1` ~ pnliprp2, data = perc)

summary(t16.1)

library(ggplot2)
ggplot(perc, aes(x = pnliprp2, y = `16:1`)) +
  geom_point(size = 4, aes(fill = treatment), shape = 21) +
  scale_fill_grey() +
  geom_line(aes(y = predict(t16.1, perc), 
                colour = "logit", linetype = "logit")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()

rlv.table.pca[,c(3:ncol(rlv.table.pca))] <- rlv.table.pca[,c(3:ncol(rlv.table.pca))] %>% mutate_if(is.character,as.numeric)

library(missMDA)
imp<-imputePCA(rlv.table.pca[,c(3:ncol(rlv.table.pca))])
imp.prop.variables<-imp$completeObs
rlv.table.pca[,c(3:ncol(rlv.table.pca))] <- imp.prop.variables

#aldh7a1 + eno3 + uggt1 + gckr + tpi1 + minpp1

test <- glm(glucose ~ eno3, data = rlv.table.pca)
summary(test)
plot(test)
car::vif(test)
