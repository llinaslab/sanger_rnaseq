df  <- data.frame(x = runif(10,min = 0,max = 10)) %>%
  tbl_df %>%
  mutate(y = x * 2 + rnorm(10, mean = 0,sd = 5)) 


df <- sample_n(bias_file, size = 1000, replace = T) %>% 
  select(exp,gc)

df <- df[which(!is.na(df$exp)),]
df <- df[which(df$gc <= 0.5),]

sampled_models <- data.frame(nrep = seq_len(1000)) %>% 
  group_by(nrep) %>% 
  do(df[df %>% nrow %>% sample.int(replace = TRUE),]) %>% 
  do(model = lm(exp ~ gc, data = .)) 

summary(unlist(lapply(sampled_models$model, function(x){coef(x)[["gc"]]})))

xs <- df$gc %>% range %>% (function(rg) seq(rg[1],rg[2],length.out = 50)) 

boot_pred <- sampled_models %>% 
  rowwise %>% 
  do(data.frame(ys = predict(.$model,list(gc = xs)))) %>%
  ungroup %>%
  cbind(xs) 

boot_pred %>% 
  mutate(run = rep(1:1000,each = 50)) %>% 
  ggplot(aes(x = xs, y = ys, group = run %>% factor)) + 
  geom_line(alpha=0.25) +
  xlab("%GC") +
  ylab("Normalized RPKM")

boot_pred %>% 
  group_by(xs) %>%
  summarize(up = quantile(ys, probs = 0.975),
            lo = quantile(ys, probs = 0.025)) %>%
  mutate(actual_data = lm(exp ~ gc, data = df) %>% predict(list(gc = xs))) %>%
  (function(d) {ggplot(df,aes(x = gc, y = exp)) + geom_point() + stat_smooth(method="lm",color="red") +
      geom_line(aes(x = xs, y = up),data = d) +
      geom_line(aes(x = xs, y = lo),data = d)
  })