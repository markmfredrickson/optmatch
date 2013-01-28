
set.seed(20130125) # so that the random values we generate will be consistent
N <- 1000
X <- data.frame(X1 = rnorm(N), 
                X2 = rnorm(N, mean = runif(N, -5, 5)), 
                X3 = as.factor(sample(letters[1:5], N, replace = T)))

mm <- model.matrix(I(rep(1, N)) ~  X1 + X2 + X1:X3, data = X)
coefs <- runif(dim(mm)[2], -2, 2)
logits <- as.vector(coefs %*% t(mm)) 
DATA <- data.frame(Z = rbinom(N, size = 1, prob = plogis(logits)), X)
model <- glm(Z ~ X1 + X2 + X1:X3, data = DATA)
predicted <- predict(model)

save(file = "setup.rda",
     N, DATA, model, predicted)
