# import numpy as np
# from scipy.optimize import minimize
# from scipy.stats import norm
# from tqdm import tqdm
# from multiprocessing import Pool




# def pairwise_overlap_ests(sketches, seqs, est_overlap, **args):
#     estimates = []
#     n = len(seqs)
#     for i in tqdm(range(n-1), desc='Estimating pairwise overlaps', leave=True):
#         for j in range(i+1,n):
#             arr = est_overlap(sketches[i], sketches[j], seqs[i], seqs[j]) 
#             estimates.append((i+1,j+1,arr,len(seqs[i]),len(seqs[j])))
#             estimates.append((j+1,i+1,arr,len(seqs[j]),len(seqs[i])))
#     return estimates

# def pairwise_overlap_ests(sketches, seqs, est_overlap):
#     estimates,args = [],[]
#     n = len(seqs)
#     for i in range(n-1):
#         for j in range(i+1,n):
#             args.append((sketches[i], sketches[j], seqs[i], seqs[j]))
#             estimates.append((i+1,j+1,None,len(seqs[i]),len(seqs[j])))
#     args = (x for x in args)
#     with Pool() as pool:
#         arrs = pool.starmap(est_overlap, args)
#     for i in range(len(arrs)):
#         x = estimates[i]
#         estimates[i] = (x[0],x[1],arrs[i],x[3],x[4])
#     return estimates

                
                
                
# def func(x, xs, log_ys):
#     mu,sigma,offset = x
#     return np.sum(
#         [(np.log(norm.pdf(x,loc=mu,scale=sigma))+np.log(offset)-log_y)**2 \
#             for x,log_y in zip(xs,log_ys)])

# def fit_normal(vals, init, bounds, sgn):
#     vals = np.sort(vals)
#     vals = vals[:len(vals)//2] if method==1 else vals[len(vals)//2:]
#     xs, ys = np.unique(vals, return_counts=True)
#     log_ys = np.log(ys/sum(ys))
#     mu,sigma,offset = minimize(func, x0=init, args=(xs,log_ys), bounds=bounds).x
#     print('mu:', mu, 'stdv:', sigma, 'offset:', offset, 'thresh:',mu+sgn*2*sigma)
#     return mu,sigma

# def get_threshold(sketches, seqs, est_overlap, **args):
#     vals = []
#     n = len(seqs)
#     for i in np.random.choice(np.arange(n), size=n//10, replace=False):
#         for j in np.random.choice(np.arange(n), size=n//20, replace=False):
#             if i==j: continue
#             _, arr = est_overlap(sketches[i], sketches[j], seqs[i], seqs[j])
#             vals.extend(arr)
#     mu,sigma = fit_normal(vals, **args)
#     thresh = mu + args['sgn']*2*sigma
#     return thresh