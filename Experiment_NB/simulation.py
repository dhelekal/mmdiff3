import tensorflow as tf
import numpy as np
import pandas as pd
import argparse
import json
import sys
import os

tfd = tf.distributions

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def createNucleosomeDraws(mod_count, hist_num, N):
    nucleosome_width = 300
    mus = tf.reshape(tf.tile(tf.linspace(nucleosome_width/2, (hist_num-1)*nucleosome_width+nucleosome_width/2, hist_num), [mod_count]), [mod_count, hist_num])
    rvs = tf.truncated_normal([N, mod_count, hist_num], mean=0.0, stddev=150.0, dtype=tf.float32)  
    location_draws = tf.add(rvs, mus)
    location_draws = location_draws - tf.reduce_min(location_draws)
    return (location_draws) ###[N x mod_count x hist_num]

def createRandomDrivers(alpha, beta, hist_num, mod_count, N): 
    nuc_probs = tfd.Beta(alpha,beta).sample([mod_count, hist_num])
    nuc_distribs = tfd.Bernoulli(probs=nuc_probs)
    samples = nuc_distribs.sample(N)
    return (samples, nuc_distribs) ###[N x mod_count x hist_num]

def createConfoundingDrivers(confounding_vector, hist_num, mod_count, N, alpha, beta):
    nuc_probs = tfd.Beta(alpha,beta).sample([hist_num]) ###[hist_num]
    nuc_distrib = tfd.Bernoulli(probs=nuc_probs) ###[hist_num]
    
    confounding_distrib = tfd.Bernoulli(probs=confounding_vector) ###[mod_count]
    
    nuc_draw = tf.reshape(nuc_distrib.sample(N), [N,1,hist_num]) ###[N x 1 x hist_num]
    confounder_sample = tf.reshape(confounding_distrib.sample(N),[N, mod_count, 1]) ###[N x mod_count x 1]
    
    samples = tf.matmul(confounder_sample, nuc_draw) ###[N x mod_count x hist_num]
    
    return (samples, confounding_distrib, nuc_distrib)
    

def createRandomGroupDrivers(confounding_vector, hist_num, mod_count, N, alpha, beta):
    nuc_probs = tfd.Beta(alpha,beta).sample([mod_count, hist_num])
    nuc_distrib = tfd.Bernoulli(probs=nuc_probs) ###[hist_num]
    nuc_draw = nuc_distrib.sample(N) ###[N x 1 x hist_num]
    
    confounding_distrib = tfd.Bernoulli(probs=confounding_vector) ###[mod_count]
    confounder_sample = tf.matrix_diag(confounding_distrib.sample(N))
    
    samples = tf.matmul(confounder_sample, nuc_draw) ###[N x mod_count x hist_num]
    
    return (samples, confounding_distrib, nuc_distrib)

def create_mod_index(mod_count, hist_num, N):
    return tf.broadcast_to(tf.range(mod_count), [N*hist_num, mod_count])  ####[N*hist_num, mod_count]

def flatten_reshape(mod_count, hist_num, N, tens): ###expect [N x mod_count x hist_num] return  [N*hist_num, mod_count]
    tens_t = tf.transpose(tens, [0, 2, 1])
    tens_r = tf.reshape(tens_t, [-1, mod_count])
    return(tens_r)
    
def draw_sample(tens, index, mask, intensity_mask=1):
    tens_f = tf.reshape(tens, [-1])
    index_f = tf.reshape(index, [-1])
    mask_f = tf.reshape(mask, [-1])
    bool_mask_f = tf.greater(mask_f, 0)
    
    index_masked = tf.boolean_mask(index_f, bool_mask_f)
    tens_masked = tf.boolean_mask(tens_f, bool_mask_f)

    return (tf.stack([tf.cast(tens_masked,dtype=tf.int32),index_masked]))
    
def create_cluster(conf_vec, mod_count, hist_num, N, alpha, beta):
    draws_fg = createNucleosomeDraws(mod_count, hist_num, N)
    driver, conf_dist, nuc_dist = createConfoundingDrivers(conf_vec, hist_num, mod_count, N, alpha, beta)

    draws = flatten_reshape(mod_count, hist_num, N, draws_fg)
    driver = flatten_reshape(mod_count, hist_num, N, driver)
    index = create_mod_index(mod_count, hist_num, N)
    
    samples = draw_sample(draws, index, driver)
    
    return(samples, nuc_dist, conf_dist)
    
def create_random_group(conf_vec, mod_count, hist_num, N, alpha, beta):
    draws_fg = createNucleosomeDraws(mod_count, hist_num, N)
    driver, conf_dist, nuc_dist = createRandomGroupDrivers(conf_vec, hist_num, mod_count, N, alpha, beta)

    draws = flatten_reshape(mod_count, hist_num, N, draws_fg)
    driver = flatten_reshape(mod_count, hist_num, N, driver)
    index = create_mod_index(mod_count, hist_num, N)
    
    samples = draw_sample(draws, index, driver)
    
    return(samples, nuc_dist, conf_dist)
    
def run_simulation( 
                    name,
                    cluster_patterns,
                    n_patterns,
                    n_marks, 
                    n_histones, 
                    n_loc_per_matrix, 
                    coverage_per_histone_bg, 
                    coverage_per_histone_fg, 
                    nuc_width, 
                    sigma_nuc,
                    alpha_pattern,
                    beta_pattern,
                    alpha_bg, 
                    beta_bg,
                    alpha_mixing_noise,
                    beta_mixing_noise, 
                    ):
    
    if not name:
        print("Error: empty name not allowed ")
        return
    
    if not os.path.exists(name+"/"):
        os.makedirs(name+"/")
        with open(name+"/patterns"+".json", "w") as text_file:
            jstr = json.dumps(cluster_patterns, indent=4, cls=NumpyEncoder)
            print(jstr, file=text_file)
        
    else:
        print("Error: directory already exits: "+name)
        return
        
    session = tf.InteractiveSession()
    
    mod_count = n_marks #Number of hist mods
    hist_num = n_histones  #Number of nucleosomes 
    N_bg = coverage_per_histone_bg #Number of Draws
    N_cf = coverage_per_histone_fg #Number of Draws 
    alpha = 0.3
    beta = 0.3
    alpha_beta = 1.
        
    for i in range(0, n_patterns):
        
        for j in range(0, len(cluster_patterns[i])):
            conf_clust = tf.constant(cluster_patterns[i][j],dtype=tf.float32)
            conf_clust = conf_clust*tfd.Beta(5.,0.8).sample([mod_count])
            clust, nc_dist, cd = create_cluster(conf_clust, mod_count, hist_num, N_cf, alpha, beta)
            with tf.name_scope('summaries'):
                tf.summary.tensor_summary("Association "+str(i)+" Cluster "+str(j), cd.probs)
            if j > 0:
                clusts = tf.concat([clust,clusts],axis=1)
            else:
                clusts = clust
                
        draws_bg = createNucleosomeDraws(mod_count, hist_num, N_bg)
        driver_bg, bg_distribs = createRandomDrivers(alpha_beta, alpha_beta, hist_num, mod_count, N_bg)
        draws_bg = flatten_reshape(mod_count, hist_num, N_bg, draws_bg)
        driver_bg = flatten_reshape(mod_count, hist_num, N_bg, driver_bg)
        index_bg = create_mod_index(mod_count, hist_num, N_bg)

        bg_samples = draw_sample(draws_bg, index_bg, driver_bg)
        
        pattern_directory = name+"/pattern_"+str(i)
        if not os.path.exists(pattern_directory):
            os.makedirs(pattern_directory)
        
        for j in range(0, n_loc_per_matrix):
            merged = tf.summary.merge_all()
            out = session.run({'bg':tf.transpose(bg_samples, [1,0]), 'fg':tf.transpose(clusts, [1,0])})
        
            bg_df_pd = pd.DataFrame({'position':out['bg'][:,0],'modification':out['bg'][:,1]})
            fg_df_pd = pd.DataFrame({'position':out['fg'][:,0],'modification':out['fg'][:,1]})

            df = pd.concat([bg_df_pd,fg_df_pd])
        
            directory = pattern_directory+"/peak_"+str(j)
            if not os.path.exists(directory):
                os.makedirs(directory)
        
            df.to_csv(directory+"/data.csv")
            
def generate_pattern(nmod, prob):
    perm = np.random.permutation(nmod)
    size_sum = 0
    clusters = []
    while size_sum < nmod:
        clust_sz = np.random.binomial(nmod-size_sum-1, prob, 1) + 1
        clust_vec = np.pad(np.ones(clust_sz), (size_sum, nmod-size_sum-clust_sz), 'constant', constant_values=(0,0))
        clust_vec = clust_vec[perm]
        clusters.append(clust_vec)
        size_sum = size_sum+clust_sz[0]
    return(clusters)

def main():    
    class C:
        pass
        
    opts = C();
    parser = argparse.ArgumentParser(description='Generates simulation data')
    parser.add_argument('config', type=str, help="Path to config file")
    parser.add_argument('--pattern', nargs=1, help="Directory containing cofounding pattern matrices", default='')
    parser.parse_args(args=sys.argv[1:], namespace=opts)
    
    if os.path.exists(opts.config):
        conf = json.load(open(opts.config,"r"))
    else:
        print("Path doesn't exist")
        exit(1)

    patterns = []
        
    if os.path.exists(opts.pattern):
        patterns = json.load(open(opts.pattern,"r"))
    else:
        for i in range(0,conf['config']['n_patterns']):
            patterns.append(generate_pattern(conf['config']['n_modifications'], conf['cluster']['prob']))

    print ("TensorFlow version: " + tf.__version__)
    tf.set_random_seed(conf['config']['random_seed'])
    
    run_simulation(
                    conf['config']['name'],
                    patterns,
                    conf['config']['n_patterns'],
                    conf['config']['n_modifications'], 
                    conf['config']['n_histones'], 
                    conf['config']['locations_per_pattern'], 
                    conf['config']['sampling_depth_background'], 
                    conf['config']['sampling_depth_foreground'],
                    conf['histone']['width'],
                    conf['histone']['sigma'],
                    conf['foreground_prior_beta']['alpha'],
                    conf['foreground_prior_beta']['beta'],
                    conf['background_prior_beta']['alpha'], 
                    conf['background_prior_beta']['beta'],
                    conf['mixing_prior_beta']['alpha'], 
                    conf['mixing_prior_beta']['beta']
                    )            

if __name__ == '__main__':
    main()