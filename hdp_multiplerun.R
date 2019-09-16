alpha=0.5
base=0.5
gamma=0.5
V=length(voc)

x_ji = word_id
t_ji=list(0)
for (j in 1:length(x_ji)){
  t_ji[[j]]=rep(1,length(x_ji[[j]]))
}
k_jt = rep(list(1),length(x_ji))  # topic for each document and table
n_jt = rep(list(0),length(x_ji)) # number of terms for each document and table
#n_jt = list(0)
for(j in 1:length(x_ji)){
  n_jt[[j]]=length(x_ji[[j]])
}

tables = rep(list(1),length(x_ji)) # list of list
n_tables = length(x_ji)

m_k = vector()
m_k[1]= length(t_ji)
n_k = vector()
n_k[1] = sum(lengths(x_ji))
#n_kv = matrix(0,nrow=1,ncol=V) # not sure of dimension
n_kv= matrix(tabulate(unlist(x_ji)),nrow=1,ncol=V,byrow=TRUE)

topics = vector() ## not sure should be vector or list
topics[1]=1

alpha_over_T_gamma = alpha / (n_tables + gamma)
Vbase = V * base
gamma_f_k_new_x_ji = gamma / V
#cur_log_base_cache = vector()
#cur_log_V_base_cache = vector()

j=1
i=1

totiter=100
n_kv_store=list(0)
#n_kv_store[[1]]=n_kv
tables_store=list(0)
#tables_store[[1]]=tables
n_jt_store=list(0)
k_jt_store=list(0)

#main function
L <- NULL
for(iter in 1:totiter){
for(j in 1:length(x_ji)){
  x_i=x_ji[[j]]
  for(i in 1:length(x_i)){
    v = x_ji[[j]][i]
    tables_j = tables[[j]]
    t_old = t_ji[[j]][i]
    if(t_old >0){
      k_old=k_jt[[j]][t_old]
      
      #decrease counters
      n_kv[k_old, v] = n_kv[k_old, v]-1
      
      n_k[k_old] = n_k[k_old]- 1
      
      n_jt[[j]][t_old] = n_jt[[j]][t_old]-1
    }
    
    if (n_jt[[j]][t_old]==0){
      #tables that all guests are gone
      tables_j=tables_j[-which(tables_j==t_old)]
      #I think we need this but why not
      tables[[j]]=tables_j
      
      m_k[k_old] = m_k[k_old]-1
      
      n_tables =  n_tables-1
      
      alpha_over_T_gamma = alpha / (n_tables + gamma)
      
      if (m_k[k_old] == 0){
        # topic (dish) that all guests are gone
        if(k_old %in% topics){
          #print('removetopic')
          #print(k_old)
          topics=topics[-(which(topics==k_old))]
        }
      }
    }
    
    #sample from the posteriors
    t_new = sampling_t(j, i, v, tables_j)
    
    # increase counters
    t_ji[[j]][i] = t_new
    
    n_jt[[j]][t_new] = n_jt[[j]][t_new]+1
    
    k_new = k_jt[[j]][t_new]
    
    n_k[k_new] = n_k[k_new]+1
    
    n_kv[k_new, v] =  n_kv[k_new, v]+1
    
    #print(c(j,i))
    
    for(q in 1:length(m_k)){
      if(m_k[q]==0){
        if(q %in% topics){
          #print('removetopic')
          #print(q)
          #topics=topics[-(which(topics==q))]
        }
      }
    }
    
  }
  
  
  for(t in tables[[j]]){
    k_old = k_jt[[j]][t]
    n_jt_temp = n_jt[[j]][t]
    
    m_k[k_old] = m_k[k_old] - 1
    
    n_k[k_old] = n_k[k_old] - n_jt_temp
    
    if(m_k[k_old] == 0){
      if(k_old %in% topics){
        #print('removetopic')
        #print(k_old)
        topics=topics[-(which(topics==k_old))]
      }
    }
    
    
    #sampling of k
    n_jtv = count_n_jtv(j, t, k_old)
    K = length(topics)
    log_p_k = rep(0, K+1) #the plus one is to accomodate a new topic if we draw one
    for (c in 1:length(topics)){
      k = topics[c]
      n_kv_k = n_kv[k,]
      n_k_k = n_k[k]
      log_p_k[c] = log_f_k_new_x_jt(n_jt,n_jtv,n_kv_k,n_k_k) + log(m_k[k])
    }    
    log_p_k[K+1] = log_f_k_new_x_jt(n_jt,n_jtv) + log(gamma)
    #print(log_p_k)
    k_new = sampling_topic(exp(log_p_k - max(log_p_k))) #the minus max is in case the log_p_k are too negative
    
    #update the counters
    k_jt[[j]][t] = k_new
    
    m_k[k_new] = m_k[k_new] + 1
    
    n_k[k_new] = n_k[k_new] + n_jt[[j]][t]
    
    for (ind in 1:length(x_ji[[j]])){
      if(t_ji[[j]][ind] != t){next}
      v_ind=x_ji[[j]][ind]
      n_kv[k_new, v_ind]  = n_kv[k_new, v_ind] + 1
    }
  }
}
  n_kv_store[[iter]]=n_kv
  tables_store[[iter]]=tables
  n_jt_store[[iter]]=n_jt
  k_jt_store[[iter]]=k_jt
  
  phi=matrix(0,dim(n_kv)[1],V)
  for(k in 1:dim(n_kv)[1]){
    phi[k,]=(n_kv[k,]+base)/(n_k[k]+Vbase)
  }
  
  l=0
  for(j in 1:length(x_ji)){
    for(i in 1:length(x_ji[[j]])){
      k_temp=k_jt[[j]][t_ji[[j]][i]]
      v_temp=rep(0,V)
      v_temp[x_ji[[j]][i]]=1
      l <- c(l,log(dmultinom(v_temp,size=1,prob=phi[k_temp,])))
    }
  }
  L[iter]=sum(l)
  print(L[iter])
}

# draw posterior topics phi_k
phi=matrix(0,dim(n_kv)[1],V)
tword=matrix(0,dim(n_kv)[1],2)
for(k in 1:dim(n_kv)[1]){
  phi[k,]=(n_kv[k,]+base)/(n_k[k]+Vbase)
  tword[k,]=(c(voc[which.max(phi[k,])],n_k[k]))
}
sort(unique(tword[,1]))



#for doc 1
print(c(voc[which.max(phi[1,])],voc[which.max(phi[43,])],voc[which.max(phi[64,])],voc[which.max(phi[3,])]))

for(j in 1:length(x_ji)){
  print(k_jt[[j]])
}

# evaluate perplexity (later)


# 2 cache functions and log_f_new_x_ji (Qiao)
cur_log_base <- function(n){
  # calculate \sum_{i=0}^{n-1} log(i + self.base)
  s=0
  for(ic in 1:n){   #0 to n-1 
    s=s+log(ic+base-1)
  }
  return(s)
}

cur_log_V_base <- function(n){
  # calculate \sum_{i=0}^{n-1} log(i + V* base)
  s=0
  for(ic in 1:n){
    s=s+log(ic+V*base-1)
  }
  return(s)
}

log_f_k_new_x_jt <- function(n_jt,n_jtv,n_kv_k=NULL,n_k_k=0){
  #when pass in, pass in (n_jt,n_jtv,n_kv[k,:],n_k[k])
  p=cur_log_V_base(n_k_k)-cur_log_V_base(n_k_k+n_jt[[j]][t])
  for(cc in 1:length(n_jtv[1,])){
    v_1=n_jtv[2,cc]
    n_1=n_jtv[1,cc]
    if(is_null(n_kv_k)==TRUE){
      n0=0
    }else{
      n0=n_kv_k[v_1]
    }
    #print(n0)
    if(n0==0){
      p=p+cur_log_base(n0+n_1)
    }else{
      p=p+cur_log_base(n0+n_1)-cur_log_base(n0)
    }
    #print(p)
  }
  return(p)
}

#Note when pass in, pass in (n_jt,n_jtv,n_kv[k,:],n_k[k])!!
#create n_jtv like:
# n_jtv= c(4,2,3)
# n_jtv_names= c(vocid[1],vocid[2],vocid[3])


# count_njtv 
# (to be used in sampling k, needs special care of njtv structure)
count_n_jtv <- function(j, t, k_old){
  #"""count n_jtv and decrease n_kv for k_old"""
  x_i = x_ji[[j]]
  t_i = t_ji[[j]]
  n_jtv_counts = vector()
  n_jtv_names =vector()
  n_jtv= rbind(n_jtv_counts,n_jtv_names) #matrix first row as counts, second as names
  for (ic in 1:length(t_i)){
    t1=t_i[ic]
    if(t1 == t){
      v = x_i[ic]
      n_kv[k_old, v] = n_kv[k_old, v] - 1
      assign('n_kv',n_kv,envir = .GlobalEnv)
      
      if(v %in% n_jtv[2,]){
        ind=which(n_jtv[2,]==v)
        n_jtv[1,ind] = n_jtv[1,ind]+1
      }else{
        n_jtv=cbind(n_jtv,c(1,v))
      }
    }
  }
  return (n_jtv)
}



# sampling t (Qiao)
sampling_t <- function(j, i, v, tables_j){
  f_k = (n_kv[, v] + base) / (n_k + Vbase)
  p_t <- NULL
  #for(t in tables_j){
  #  p_t[t]=n_jt[[j]][t] * f_k[k_jt[[j]][t]]
  #}
  for(int in 1:length(tables_j)){
    p_t[int]=n_jt[[j]][tables_j[int]] * f_k[k_jt[[j]][tables_j[int]]]
  }

  p_x_ji = sum(m_k * f_k) + gamma_f_k_new_x_ji
  p_t=c(p_t, p_x_ji * alpha_over_T_gamma)
  
  p_t=p_t/sum(p_t)
  drawing = which.max(rmultinom(1, 1, p_t))
  
  if(drawing <= length(tables_j)){
    return(tables_j[drawing])
  }else{
    return(new_table(j, i, f_k))
  }
}


# new table (Qiao)
new_table <- function(j,i,f_k){
  # search a spare table ID
  T_j = length(n_jt[[j]])
  t_new=1
  for(t_new in 1:T_j){
    if (t_new %in% tables[[j]]==FALSE){
      break}
  }
  
  #if(t_new==T_j & t_new %in% test){
  if(all(seq(1,T_j,by=1) %in% tables[[j]])){
    # new table ID (no spare)
    t_new = T_j+1
    n_jt[[j]]=c(n_jt[[j]],0)
    assign('n_jt',n_jt,envir = .GlobalEnv)
    k_jt[[j]]=c(k_jt[[j]],0)
    assign('k_jt',k_jt,envir = .GlobalEnv)
  }
  
  tables[[j]]=c(tables[[j]],t_new)
  assign('tables',tables,envir = .GlobalEnv)
  
  n_tables = n_tables+ 1
  assign('n_tables',n_tables,envir = .GlobalEnv)
  
  alpha_over_T_gamma = alpha / (n_tables + gamma)
  assign('alpha_over_T_gamma',alpha_over_T_gamma,envir = .GlobalEnv)
  
  # sampling of k for new topic(= dish of new table)
  p_k <- vector()
  for(k in topics){
    p_k[k] = m_k[k] * f_k[k]
  }
  p_k=c(p_k,gamma_f_k_new_x_ji)
  p_k=na.omit(p_k)##########################################################
  
  k_new = sampling_topic(p_k)
  
  k_jt[[j]][t_new] = k_new
  assign('k_jt',k_jt,envir = .GlobalEnv)
  
  m_k[k_new] = m_k[k_new]+ 1
  assign('m_k',m_k, envir = .GlobalEnv)
  
  return (t_new)
  
}



# sampling topic (Jingjing)
sampling_topic <- function(p_k){
  drawing = which.max(rmultinom(1, 1, prob = p_k/sum(p_k)))
  if (drawing <= length(topics)){
    k_new = topics[drawing]
  }else{
    K = length(m_k)
    for(k_new in 1:K){
      if(k_new %in% topics == FALSE){break
      }}
    
    if(all(seq(1,K,by=1) %in% topics)){
      k_new = K+1
      n_k = c(n_k, 0)
      assign('n_k',n_k,envir = .GlobalEnv)
      
      m_k = c(m_k, 0)
      assign('m_k',m_k,envir = .GlobalEnv)
      
      n_kv = rbind(n_kv, rep(0,V))
      assign('n_kv',n_kv,envir = .GlobalEnv)
    }
    topics = c(topics, k_new)
    assign('topics',topics,envir = .GlobalEnv)
  }
  return(k_new)
}






for(j in 1:length(x_ji)){
for(it in 1:length(x_ji[[j]])){
  t_it=t_ji[[j]][it]
  if(k_jt[[j]][t_it]==6){
    print(c(j,it))
  }
}
}


topics=c(1,3,5,7,8)
m_k=c(1,0,1,1,1)
for(q in 1:length(m_k)){
  if(m_k[q]==0){
    if(q %in% topics){
      topics=topics[-(which(topics==q))]
    }
  }
}

save.image("~/Desktop/StatModeling2/Project-HDP/testrun2-success10iter.RData")
