void load_convolutional_weights(layer l, FILE *fp)
{
    if(l.binary){
        //load_convolutional_weights_binary(l, fp);
        //return;
    }
    if(l.numload) l.n = l.numload;
    int num = l.c/l.groups*l.n*l.size*l.size;
    fread(l.biases, sizeof(float), l.n, fp);
    if (l.batch_normalize && (!l.dontloadscales)){
        fread(l.scales, sizeof(float), l.n, fp);
        fread(l.rolling_mean, sizeof(float), l.n, fp);
        fread(l.rolling_variance, sizeof(float), l.n, fp);
        if(0){
            int i;
            for(i = 0; i < l.n; ++i){
                printf("%g, ", l.rolling_mean[i]);
            }
            printf("\n");
            for(i = 0; i < l.n; ++i){
                printf("%g, ", l.rolling_variance[i]);
            }
            printf("\n");
        }
        if(0){
            fill_cpu(l.n, 0, l.rolling_mean, 1);
            fill_cpu(l.n, 0, l.rolling_variance, 1);
        }
        if(0){
            int i;
            for(i = 0; i < l.n; ++i){
                printf("%g, ", l.rolling_mean[i]);
            }
            printf("\n");
            for(i = 0; i < l.n; ++i){
                printf("%g, ", l.rolling_variance[i]);
            }
            printf("\n");
        }
    }
    fread(l.weights, sizeof(float), num, fp);
    //if(l.c == 3) scal_cpu(num, 1./256, l.weights, 1);
    if (l.flipped) {
        transpose_matrix(l.weights, l.c*l.size*l.size, l.n);
    }
    //if (l.binary) binarize_weights(l.weights, l.n, l.c*l.size*l.size, l.weights);
#ifdef GPU
    if(gpu_index >= 0){
        push_convolutional_layer(l);
    }
#endif

///add start
    int i,j;
    FILE *fp_w = fopen("weights.bin", "ab+");
    if(!fp_w) file_error("weights.bin");
    FILE *fp_bias = fopen("bias.bin", "ab+");
    if(!fp_bias) file_error("bias.bin");

    float *weight_buffer = (float *)calloc(num, sizeof(float));
    float *alpha_buffer = (float *)calloc(l.n, sizeof(float));
    float *bias_buffer = (float *)calloc(l.n, sizeof(float));

    if(l.batch_normalize && (!l.dontloadscales))
    {
        for(i = 0;i < l.n; i++)
        {
            float tmp = l.scales[i]/(sqrt(l.rolling_variance[i]) + .000001f);
            alpha_buffer[i] = tmp;
            bias_buffer[i] = l.biases[i] - l.rolling_mean[i]*tmp;
        }
    }
    else
    {
        for(i = 0;i < l.n; i++)
        {
            alpha_buffer[i] = 1;
            bias_buffer[i] = l.biases[i];
        }
    }

    int cnt = 0;
    for(j = 0;j < l.n; j++)
        for(i = 0;i < num/l.n; i++)
        {
            weight_buffer[cnt] = l.weights[cnt]*alpha_buffer[j];
            cnt++;
        }

    fwrite(weight_buffer, sizeof(float), num, fp_w);
    fwrite(bias_buffer, sizeof(float), l.n, fp_bias);

    fclose(fp_w);
    fclose(fp_bias);

    free(weight_buffer);
    free(alpha_buffer);
    free(bias_buffer);    
///add end 

}
