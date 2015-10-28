#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FFT_SIZE 512
#define TWICEPI (6.2831852f)
#define	PI      (3.1415926f)

typedef struct{
    float real;
    float imag;
}complex;

static complex multicomplex(complex b1, complex b2) /* multiplication of complex */
{
    complex b3;
    b3.real=b1.real*b2.real-b1.imag*b2.imag;
    b3.imag=b1.real*b2.imag+b1.imag*b2.real;
    return(b3);
}

static int mylog2(int N) /* Max(N) = 4098 */
{
    int k=0;

    if (N>>12) { k+=12; N>>=12; }
    if (N>>8) { k+=8; N>>=8; }
    if (N>>4) { k+=4; N>>=4; }
    if (N>>2) { k+=2; N>>=2; }
    if (N>>1) { k+=1; N>>=1; }
    return k ;
}

static void BitReverse(complex *xin, int N)
{
    int LH, i, j, k;
    complex tmp;

    LH=N/2;
    j = N/2;

    for( i = 1; i <= (N -2); i++)
    {
        if(i < j)
        {
            tmp = xin[j];
            xin[j] = xin[i];
            xin[i] = tmp;
        }
        k = LH;
        while(j >= k)
        {
            j = j-k;
            k = k/2;
        }
        j = j + k;
    }
}

static void DFT_2(complex *b1, complex *b2)
{
    complex tmp;
    tmp = *b1;

    (*b1).real = (*b1).real + (*b2).real;
    (*b1).imag = (*b1).imag + (*b2).imag;

    (*b2).real = tmp.real - (*b2).real;
    (*b2).imag = tmp.imag - (*b2).imag;
}

static void DFT_4(complex* b0, complex* b1, complex* b2, complex* b3)
{
    /*variables locales*/
    complex temp[4];

    /*calcul x1*/
    temp[0].real=(*b0).real+(*b1).real;
    temp[0].imag=(*b0).imag+(*b1).imag;

    /*calcul x2*/
    temp[1].real=(*b0).real-(*b1).real;
    temp[1].imag=(*b0).imag-(*b1).imag;

    /*calcul x3*/
    temp[2].real=(*b2).real+(*b3).real;
    temp[2].imag=(*b2).imag+(*b3).imag;

    /*calcul x4 + multiplication with -j*/
    temp[3].imag=(*b3).real-(*b2).real;
    temp[3].real=(*b2).imag-(*b3).imag;

    /*the last stage*/
    (*b0).real=temp[0].real+temp[2].real;
    (*b0).imag=temp[0].imag+temp[2].imag;

    (*b1).real=temp[1].real+temp[3].real;
    (*b1).imag=temp[1].imag+temp[3].imag;

    (*b2).real=temp[0].real-temp[2].real;
    (*b2).imag=temp[0].imag-temp[2].imag;

    (*b3).real=temp[1].real-temp[3].real;
    (*b3).imag=temp[1].imag-temp[3].imag;
}



static void FFT_R4(complex *xin, int N, int m)
{
    int i, L, j;
    double ps1, ps2, ps3;
    int le,B;
    complex w[4];

    for( L = 1; L <= m; L++)
    {
        le = pow(4 ,L);
        B = le/4; /*the distance of buttefly*/

        for(j = 0; j <= B-1 ; j++)
        {
            // ps0 = (TWICEPI/N) * 0 * j;
            // w[0].real = cos(ps0);
            // w[0].imag = -sin(ps0);

            ps1 = ((TWICEPI)/le)*2*j;
            w[1].real = cos(ps1);
            w[1].imag = -sin(ps1);

            ps2 = (TWICEPI/le)*j;
            w[2].real = cos(ps2);
            w[2].imag = -sin(ps2);

            ps3 = (TWICEPI/le)*3*j;
            w[3].real = cos(ps3);
            w[3].imag = -sin(ps3);

            for(i = j; i <= N-1; i = i + le) /* controle those same butteflies*/
            {
                /* multiple with W */
                // xin[i] = multicomplex(xin[i], w[0]);
                xin[i + B] = multicomplex(xin[i + B], w[1]);
                xin[i + 2*B] = multicomplex(xin[i + 2*B], w[2]);
                xin[i + 3*B] = multicomplex(xin[i + 3*B], w[3]);
                /* DFT-4 */
                DFT_4(xin + i, xin + i + B, xin + i + 2*B, xin + i + 3*B);
            }
        }
        /*
           printf("*****N°%d **********\n", L);
           for(i=0;i<N;i++)
           {
           printf("%.8f\t\t",xin[i].real);
           printf("%.8f\n",xin[i].imag);
           }
           */
    }

}//end of FFT_R4

static void FFT_L2(complex *xin, int N)
{ /* For the last stage 2 DFT*/
    int j, B;
    double p, ps ;
    complex w;

    B = N/2;
    for(j = 0; j <= B - 1; j++)
    {
        ps = (TWICEPI/N)*j;
        w.real = cos(ps);
        w.imag = -sin(ps);

        /* multiple avec W */
        xin[j+ B] = multicomplex(xin[j + B], w);
        DFT_2(xin + j ,xin + j + B);
    }

}//end of FFT_L2


/**************** FFT complex ***************/
void RunFFT(complex *xin, int N)
{
    int m, i;
    BitReverse(xin, N);
    m = mylog2(N);
    if( (m%2) == 0 )
    {
        /*All the stages are radix 4*/
        FFT_R4(xin, N, m/2);
    }
    else
    {
        /*the last stage is radix 2*/
        FFT_R4(xin, N, m/2);
        FFT_L2(xin, N);
    }
}

void RunIFFT(complex *xin,int N)
{ /* inverse FFT */
    int i;

    for(i=0; i < N + 1 ; i++)
    {
        xin[i].imag = -xin[i].imag;
    }

    RunFFT(xin,N);

    for(i = 0; i < N + 1 ; i++)
    {
        xin[i].real = xin[i].real/N;
        xin[i].imag = -xin[i].imag/N;
    }
}

/************** FFT real ****************/
void RunFFTR(complex *xin, int N)
{
    int i;
    double ps;
    complex *Realin;
    complex Realtmp1;
    complex Realtmp2;
    complex w;

    Realin = (complex *)malloc((N/2)*sizeof(complex));

    /*** X(n)= A(2n) +j*A(2n+1) ***/
    for(i = 0 ; i < N/2 ; i++)
    {
        Realin[i].real = xin[2*i].real;
        Realin[i].imag = xin[2*i + 1].real;
    }

    RunFFT(Realin, N/2);

    for(i = 0; i < N/2 ; i++)
    {
        /***** factor w *****/
        ps = (TWICEPI/N)*i;
        w.real = cos(ps);
        w.imag = -sin(ps);

        /***** conjugue *****/
        Realtmp1.real = Realin[(N/2) - i].real;
        Realtmp1.imag = -Realin[(N/2) - i].imag;

        if(i == 0)
        {
            Realtmp1.real = Realin[0].real;
            Realtmp1.imag = -Realin[0].imag;
        }
        /***** part 2 *****/
        Realtmp2.real = Realin[i].imag - Realtmp1.imag;
        Realtmp2.imag = Realtmp1.real - Realin[i].real;
        if(i > 0)
            Realtmp2 = multicomplex(Realtmp2, w);

        /**** part 1 ****/
        Realtmp1.real = Realin[i].real + Realtmp1.real;
        Realtmp1.imag = Realin[i].imag + Realtmp1.imag;

        xin[i].real = (Realtmp1.real + Realtmp2.real)/2;
        xin[i].imag = (Realtmp1.imag + Realtmp2.imag)/2;

        xin[N/2 + i].real = (Realtmp1.real - Realtmp2.real)/2;
        xin[N/2 + i].imag = (Realtmp1.imag - Realtmp2.imag)/2;
    }
}

void RunIFFTR(complex *xin, int N)
{

    int i;
    double ps;
    complex *Realin;
    complex Realtmp1;
    complex Realtmp2;
    complex w;

    Realin = (complex *)malloc((N/2)*sizeof(complex));

    /*** X(n)= A(2n) +j*A(2n+1) ***/
    for(i = 0 ; i < N/2 ; i++)
    {
        Realin[i].real = xin[2*i].real;
        Realin[i].imag = xin[2*i + 1].real;
    }

    RunIFFT(Realin, N/2);

    for(i = 0; i < N/2 ; i++)
    {
        /***** factor w *****/
        ps = (TWICEPI/N)*i;
        w.real = cos(ps);
        w.imag = -sin(ps);

        /***** conjugue *****/
        Realtmp1.real = Realin[(N/2) - i].real;
        Realtmp1.imag = -Realin[(N/2) - i].imag;

        if(i == 0)
        {
            Realtmp1.real = Realin[0].real;
            Realtmp1.imag = -Realin[0].imag;
        }
        /***** part 2 *****/
        Realtmp2.real = Realin[i].imag - Realtmp1.imag;
        Realtmp2.imag = Realtmp1.real - Realin[i].real;
        if(i > 0)
            Realtmp2 = multicomplex(Realtmp2, w);

        /**** part 1 ****/
        Realtmp1.real = Realin[i].real + Realtmp1.real;
        Realtmp1.imag = Realin[i].imag + Realtmp1.imag;

        xin[i].real = (Realtmp1.real + Realtmp2.real)/2;
        xin[i].imag = (Realtmp1.imag + Realtmp2.imag)/2;

        xin[N/2 + i].real = (Realtmp1.real - Realtmp2.real)/2;
        xin[N/2 + i].imag = (Realtmp1.imag - Realtmp2.imag)/2;

    }
}

void FFT_PCM(float *xout, float *xin, short fft_size)
{
    int i;
    complex x[1024];
    for(i=0; i<fft_size; i++)
    {
        x[i].real = (float)xin[i];
        x[i].imag = 0;
    }
    RunFFTR(x, fft_size);
    for(i=0; i<fft_size/2; i++)
        xout[i] = sqrt(x[i].real*x[i].real + x[i].imag*x[i].imag);
}

/* Supporting routine for MFCC */
#define MfccRound(x) ((int)((x)+0.5))
#define BANK_NUM 26
typedef struct {
    int m_lowX;
    int m_centerX;
    int m_highX;
    float m_sumWeight;
} WfMelFB; /* mel filter bank for noise reduction */
WfMelFB m_MelFB[BANK_NUM];
float m_MelWeight[FFT_SIZE/2+1];
float m_dctMatrix[(12+1)*BANK_NUM];

//初始化 DCT 矩阵，查表方式调用 DCT ，可以优化
int MfccInitDCTMatrix (float *dctMatrix, int ceporder, int numChannels)
{
    int i, j;
    for (i = 0; i <= ceporder; i++)      //12+1
    {
        for (j = 0; j < numChannels; j++)  //BANK_NUM
        {
            dctMatrix[i * numChannels + j] = (float)cos(PI * (float) i / (float) numChannels * ((float) j + 0.5));
            if(i==0) 
                dctMatrix[i * numChannels + j] *= (float)sqrt(1/(float)numChannels);
            else     
                dctMatrix[i * numChannels + j] *= (float)sqrt(2/(float)numChannels);
        }
    }
    return 1;
}

//初始化Mel滤波器组，程序启动时运行，将滤波器参数存入全局数组，多线程程序需要优化
void MfccInitMelFilterBanks (float startingFrequency, float samplingRate, int fftLength, int numChannels)
{
    int i, k;
    float* freq=(float*)malloc(numChannels*4+8);
    int * bin=(int *)malloc(numChannels*4+8);
    float start_mel, fs_per_2_mel;

    //m_MelWeight = (float*)malloc(fftLength*2+4);
    /* Constants for calculation */
    freq[0] = startingFrequency;
    start_mel = (float)(2595.0 * log10 (1.0 + startingFrequency / 700.0));
    bin[0] = MfccRound(fftLength * freq[0] / samplingRate);
    freq[numChannels+1] = (float)(samplingRate / 2.0);
    fs_per_2_mel = (float)(2595.0 * log10 (1.0 + (samplingRate / 2.0) / 700.0));
    bin[numChannels+1] = MfccRound(fftLength * freq[numChannels+1] / samplingRate);

    /* Calculating mel-scaled frequency and the corresponding FFT-bin */
    /* number for the lower edge of the band                          */
    for (k = 1; k <= numChannels; k++) 
    {
        freq[k] = (float)(700 * (pow (10, (start_mel + (float) k / (numChannels + 1) * (fs_per_2_mel - start_mel)) / 2595.0) - 1.0));
        bin[k] = MfccRound(fftLength * freq[k] / samplingRate);
    }

    /* This part is never used to compute MFCC coefficients */
    /* but initialized for completeness                     */


    for(i = 0; i<bin[0]; i++){
        m_MelWeight[i]=0;
    }
    m_MelWeight[fftLength/2]=1;

    /* Initialize low, center, high indices to FFT-bin */
    for (k = 0; k <= numChannels; k++) 
    {
        if(k<numChannels)
        {
            m_MelFB[k].m_lowX=bin[k];
            m_MelFB[k].m_centerX=bin[k+1];
            m_MelFB[k].m_highX=bin[k+2];
        }
        for(i = bin[k]; i<bin[k+1]; i++)
        {
            m_MelWeight[i]=(i-bin[k]+1)/(float)(bin[k+1]-bin[k]+1);
            //printf("%f,",m_MelWeight[i]);
        }
        //printf("\n");
    }


    free(freq);
    free(bin);


    return;
}
//输出mel滤波后的频域能量，滤波加降维
void MfccMelFilterBank (float *sigFFT, int numChannels, float* output, int normalize)
{
    float sum, wsum;
    int i, k;
    WfMelFB   *melFB;

    for (k=0;k<numChannels;k++)
    {
        melFB = m_MelFB+k;
        sum = sigFFT[melFB->m_centerX];
        wsum=1;
        for (i = melFB->m_lowX; i < melFB->m_centerX; i++)
        {
            sum += m_MelWeight[i] * sigFFT[i];
            wsum += m_MelWeight[i];
        }
        for (i = melFB->m_centerX+1; i <= melFB->m_highX; i++)
        {
            sum += (1 - m_MelWeight[i-1]) * sigFFT[i];
            wsum += (1 - m_MelWeight[i-1]);
        }
        output[k] = sum;
        if(normalize) 
        {
            output[k] /= wsum;
        }
    }
    return;
}

//PCM预处理,会将short型的pcm转换成float
int preemphasize(float *sample, short *pcm_in, int sampleN)
{
    /* Setting emphFac=0 turns off preemphasis. */
    int	i;
    float emphFac = (float)0.97;

    for (i = sampleN-1; i > 0; i--) {
        sample[i] = (float)pcm_in[i] - emphFac * pcm_in[i-1];
    }
    sample[0] = (float)(1.0 - emphFac) * sample[0];

    return(1);
}
//计算 DCT 需要优化
void MfccDCT (float *mel_cep, float *x)
{
    int i, j;
    for (i = 0; i <= 12; i++) 
    {
        mel_cep[i] = 0.0;
        for (j = 0; j < BANK_NUM; j++)
        {
            mel_cep[i] += x[j] * m_dctMatrix[i * BANK_NUM + j];
        }
    }	
    return;
}
//mfcc特征提取接口函数
int mfcc_extract(float *feat_out, short *pcm_in, short fft_size)
{
    int i;
    float fft_data[FFT_SIZE] = {0.0};
    float pcm_out[FFT_SIZE] = {0.0};
    float mel_bank_out[26] = {0.0};
    float energy = 0.0;
    preemphasize(pcm_out, pcm_in, FFT_SIZE);
#if 0
    for(i=0; i<512; i++)
        printf("%f,", pcm_out[i]);
#endif
    FFT_PCM(fft_data, pcm_out, 512);
    for(i=0; i<FFT_SIZE/2+1; i++)
    {
        fft_data[i] = (1.0/FFT_SIZE) * (fft_data[i]) * (fft_data[i]);
        energy += fft_data[i];
        //printf("%f,",fft_data[i]);
    }
    //printf("\n");
    if(energy == 0)
        energy = 0.0001;
    MfccMelFilterBank(fft_data, BANK_NUM, mel_bank_out, 1);
    for(i=0; i<BANK_NUM; i++)
    {
        if(mel_bank_out[i] == 0)
            mel_bank_out[i] = 0.0001;
        mel_bank_out[i] = log(mel_bank_out[i]);
        //printf("%f,",mel_bank_out[i]);
    }
    //printf("\n");
    MfccDCT(feat_out, mel_bank_out);

    return 1;
}

#if 1
int main()
{
    short data[512] = {0};
    float feat_out[12] = {0};
    short i;

    MfccInitDCTMatrix (m_dctMatrix, 12, BANK_NUM);
    MfccInitMelFilterBanks(0, 16000, FFT_SIZE, BANK_NUM);
#if 0
    for(i=0; i<BANK_NUM; i++)
        printf("%d,%d,%d\n",m_MelFB[i].m_lowX,m_MelFB[i].m_centerX,m_MelFB[i].m_highX);
#endif
    for(i=0; i<512; i++)
        data[i] = i;
    mfcc_extract(feat_out, data, 512);
#if 1
    printf("mel_cep:\n");
    for(i=0; i<13; i++)
        printf("%f,", feat_out[i]);
#endif
}
#endif
