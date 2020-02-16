#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define Q31_1_BASE (1LL<<30)
#define Q31_1_MAX  ((1<<30)-1)
#define Q31_1_MIN  (-1<<30)
#define SIZE (128)
#define PI (3.14159265358979323846)

int32_t flt2fixd(double x)
{
	if (x >= 1)
		return Q31_1_MAX;
	else if (x < -1)
		return Q31_1_MIN;

	int32_t res = x * (double)Q31_1_BASE;
	return res;
}

float fixd2flt(int32_t x)
{
	float res = (float)(x) / ((float)Q31_1_BASE);
	return res;
	return res;
}

typedef struct
{
	uint8_t ChunkID[4];
	uint32_t ChunkSize;
	uint8_t Format[4];
	uint8_t Subchunk1ID[4];
	uint32_t Subchunk1Size;
	uint16_t AudioFormat;
	uint16_t NumChannels;
	uint32_t SampleRate;
	uint32_t ByteRate;
	uint16_t BlockAlign;
	uint16_t BitsPreSample;

	uint8_t ID[4];
	uint32_t size;
} wav_header;

void wav_header_init(wav_header* header, int32_t sample_lenth) {
	header->ChunkID[0] = 'R';
	header->ChunkID[1] = 'I';
	header->ChunkID[2] = 'F';
	header->ChunkID[3] = 'F';
	header->Format[0] = 'W';
	header->Format[1] = 'A';
	header->Format[2] = 'V';
	header->Format[3] = 'E';
	header->Subchunk1ID[0] = 'f';
	header->Subchunk1ID[1] = 'm';
	header->Subchunk1ID[2] = 't';
	header->Subchunk1ID[3] = ' ';

	header->AudioFormat = 1; //pcm
	header->NumChannels = 2;
	header->SampleRate = 48000;
	header->BitsPreSample = 16;
	header->ByteRate = (header->SampleRate * header->BitsPreSample * header->NumChannels) / 8; //!
	header->BlockAlign = header->NumChannels * header->BitsPreSample / 8;
	header->Subchunk1Size = 16;
	header->size = sample_lenth * header->NumChannels * header->BitsPreSample / 8;
	header->ID[0] = 'd';
	header->ID[1] = 'a';
	header->ID[2] = 't';
	header->ID[3] = 'a';
	header->ChunkSize = 4 + (8 + header->Subchunk1Size) + (8 + header->size);
}


int32_t sine_wave(float  amplitude, float fq, int t, int32_t fs) {
	int32_t sample_sig;
	float omeg = 2 * PI * fq;
	sample_sig = flt2fixd(amplitude * sin(omeg * (float)t / fs));

	return (sample_sig >> 15);
}

int32_t noise(float  amplitude, int t, int32_t fs) {
	int32_t sample_sig;
	float omeg = 2 * PI;
	sample_sig = flt2fixd(amplitude * sin(rand() * omeg * (float)t / fs));

	return (sample_sig >> 15);
}

int16_t sweep_wave(float amplitude, float start_fq, float end_fq, int32_t sample_lenth, int t, int32_t fs)
{
	int32_t sample_sig;
	float omeg1 = 2 * PI * start_fq;
	float omeg2 = 2 * PI * end_fq;
	float tmp1 = log(omeg2 / omeg1);

	float n;
	float tmp2;
	n = (float)t / sample_lenth;
	tmp2 = exp(n * tmp1) - 1.0;
	sample_sig = flt2fixd(amplitude * sin(omeg1 * sample_lenth * tmp2 / (fs * tmp1)));

	return sample_sig >> 15;
}

int64_t acc = 0;

int32_t IIR(int32_t* buffer, int32_t* coeffs, int16_t sample)
{
	buffer[0] = buffer[1];
	buffer[1] = buffer[2];
	buffer[2] = (int32_t)sample;
	buffer[3] = buffer[4];
	buffer[4] = buffer[5];

	acc += (int64_t)buffer[2] * coeffs[0] + (int64_t)buffer[1] * coeffs[1] + (int64_t)buffer[0] * coeffs[2] - (int64_t)buffer[4] * coeffs[3] - (int64_t)buffer[3] * coeffs[4];

	buffer[5] = (acc >> 30);
	acc = acc & 0x3fffffff;
	return buffer[5];
}

void main() {
	int32_t fs = 48000;
	float freq = 30;
	float end_freq = 20000;
	float amplitude = 0.8;
	int32_t sample_lenth = 0;

	int16_t output;
	int16_t outputR;

	uint32_t second = 0;
	printf("Enter duration of record:");
	scanf_s("%d", &second);
	sample_lenth = second * fs;
	printf("Amount of samples=%d\n", sample_lenth);

	FILE* file_out;
	wav_header header;
	wav_header_init(&header, sample_lenth);

	fopen_s(&file_out, "test_signal.wav", "wb");
	fwrite(&header, sizeof(header), 1, file_out);

	int16_t buffer[2*SIZE];
	int32_t sample_buffer_L[6] = { 0, 0, 0, 0, 0, 0 };
	int32_t sample_buffer_R[6] = { 0, 0, 0, 0, 0, 0 };

	int16_t signal;
	int16_t signalR;

	double a0, a1, a2, b1, b2;
	int32_t coeffs[5] = { 0, 0, 0, 0, 0 };

	a0 = 0.11205483175086794;
	a1 = 0.2241096635017359;
	a2 = 0.11205483175086794;
	b1 = -0.8559866467934273;
	b2 = 0.3042059737968993;

	coeffs[0] = flt2fixd(a0);
	coeffs[1] = flt2fixd(a1);
	coeffs[2] = flt2fixd(a2);
	coeffs[3] = flt2fixd(b1);
	coeffs[4] = flt2fixd(b2);
	int j;
	for (int t = 0; t < sample_lenth;) {

		for (int n = 0; n < SIZE; n++, t++) {

			output = sweep_wave(amplitude, freq, end_freq, sample_lenth, t, fs);
			
			signal = IIR(sample_buffer_L, coeffs, output);
			signalR = IIR(sample_buffer_R, coeffs, output);

			buffer[2*n] = signal;
			buffer[2 * n + 1] =signalR;
			j = n;
		}

		fwrite(buffer, 2*(j + 1)*sizeof(int16_t), 1, file_out);
	}

	fclose(file_out);
	printf("\n WAV File was generated.\n");
	system("pause");
}