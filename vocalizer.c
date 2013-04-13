#include <math.h>
#include <fftw3.h>
#include <string.h>
#include <stdlib.h>
#include <err.h>
#include <stdio.h>
#include <audiofile.h>
#include <unistd.h>

#define bound(x,y,z) (((y)<(x))?(x):((y)>(z))?(z):(y))

typedef float sample_t;
#define AF_SAMPFMT_sample_t AF_SAMPFMT_FLOAT

double do_fftcompare(int n, double *a, double *b)
{
	fftw_complex c_a[n/2+1];
	fftw_complex c_b[n/2+1];
	fftw_plan p_a = fftw_plan_dft_r2c_1d(n, a, c_a, FFTW_ESTIMATE);
	fftw_plan p_b = fftw_plan_dft_r2c_1d(n, b, c_b, FFTW_ESTIMATE);
	fftw_execute(p_a);
	fftw_execute(p_b);
	double r = 0;
	int i;
	for(i = 0; i < n/2+1; ++i)
	{
		double l_a = c_a[i][0] * c_a[i][0] + c_a[i][1] * c_a[i][1];
		double l_b = c_b[i][0] * c_b[i][0] + c_b[i][1] * c_b[i][1];
		r += (l_a - l_b) * (l_a - l_b);
	}
	return r;
}

double do_fftcompare_samples(int n, int chans, sample_t *a, sample_t *b)
{
	double sa[n], sb[n];
	int i, c;
	double r = 0;
	for(c = 0; c < chans; ++c)
	{
		for(i = 0; i < n; ++i)
		{
			sa[i] = a[chans * i + c];
			sb[i] = b[chans * i + c];
		}
		r += do_fftcompare(n, sa, sb);
	}
	return r;
}

double do_directcompare(int n, double *a, double *b)
{
	double r = 0;
	int i;
	for(i = 0; i < n; ++i)
		r += (a[i] - b[i]) * (a[i] - b[i]);
	return r;
}

double do_directcompare_samples(int n, int chans, sample_t *a, sample_t *b)
{
	double sa[n], sb[n];
	int i, c;
	double r = 0;
	for(c = 0; c < chans; ++c)
	{
		for(i = 0; i < n; ++i)
		{
			sa[i] = a[chans * i + c];
			sb[i] = b[chans * i + c];
		}
		r += do_directcompare(n, sa, sb);
	}
	return r;
}

void do_fftdiff(int n, double rate, double overdrive, double *in, double *kar, double *out)
{
	(void) rate;
	fftw_complex c_in[n/2+1];
	fftw_complex c_kar[n/2+1];
	fftw_complex c_out[n/2+1];
	fftw_plan p_in = fftw_plan_dft_r2c_1d(n, in, c_in, FFTW_ESTIMATE);
	fftw_plan p_kar = fftw_plan_dft_r2c_1d(n, kar, c_kar, FFTW_ESTIMATE);
	fftw_plan p_out = fftw_plan_dft_c2r_1d(n, c_out, out, FFTW_ESTIMATE);
	fftw_execute(p_in);
	fftw_execute(p_kar);
	int i;
	for(i = 0; i < n/2+1; ++i)
	{
		double l_in = c_in[i][0] * c_in[i][0] + c_in[i][1] * c_in[i][1];
		double l_kar = c_kar[i][0] * c_kar[i][0] + c_kar[i][1] * c_kar[i][1];
		l_kar *= overdrive;
		if(l_in <= l_kar)
		{
			c_out[i][0] = c_out[i][1] = 0;
		}
		else
		{
			double f = 1 - sqrt(l_kar / l_in);
			c_out[i][0] = c_in[i][0] * f / n;
			c_out[i][1] = c_in[i][1] * f / n;
		}
	}
	fftw_execute(p_out);
}

double fft_window(double in)
{
	return in;
}

double out_window(double in)
{
	return 3 * in * in - 2 * in * in * in;
}

void copy_windowed(int win, int lap, int chan, const sample_t *in, double *out)
{
	// out[i] = in[i+shift]
	int i;
	for(i = 0; i < lap; ++i)
		out[i] = in[chan * i] * fft_window((i+0.5) / (double) (lap));
	for(i = lap; i < lap + win; ++i)
		out[i] = in[chan * i];
	for(i = lap + win; i < lap + win + lap; ++i)
		out[i] = in[chan * i] * fft_window((lap + win + lap - i-0.5) / (double) (lap));
}

void add_windowed(int win, int lap, int chan, const double *in, sample_t *out)
{
	int i;
	for(i = 0; i < lap; ++i)
		out[chan * i] += in[i] * (out_window((i+0.5) / (double) (2 * lap)) / fft_window((i+0.5) / (double) (lap)));
	for(i = lap; i < lap + lap; ++i)
		out[chan * i] += in[i] * (out_window((i+0.5) / (double) (2 * lap)));
	for(i = 2 * lap; i < lap + win - lap; ++i)
		out[chan * i] += in[i];
	for(i = lap + win - lap; i < lap + win; ++i)
		out[chan * i] += in[i] * ((1 - out_window((i+0.5 - lap - win + lap) / (double) (2 * lap))));
	for(i = lap + win; i < lap + win + lap; ++i)
		out[chan * i] += in[i] * ((1 - out_window((i+0.5 - lap - win + lap) / (double) (2 * lap))) / fft_window((lap + win + lap - i-0.5) / (double) (lap)));
}

void do_filter(int win, int lap, int chan, double rate, double overdrive, const sample_t *in, const sample_t *kar, sample_t *out)
{
	int c;
	double thisin[lap+win+lap];
	double thiskar[lap+win+lap];
	double thisout[lap+win+lap];

	for(c = 0; c < chan; ++c)
	{
		copy_windowed(win, lap, chan, in + c, thisin);
		copy_windowed(win, lap, chan, kar + c, thiskar);
		do_fftdiff(lap+win+lap, rate, overdrive, thisin, thiskar, thisout);
		add_windowed(win, lap, chan, thisout, out + c);
	}
}

#define USAGE \
	"Usage: %s -i infile.wav -k karaokefile.wav -o outfile.wav\n" \
	"  [-w windowsize] [-l overlap]\n" \
	"  [-d delta_min] [-D delta_max] [-s search_min] [-S search_max]\n" \
	"  [-f overdrive]\n" \
	""

typedef struct
{
	const char *infile;
	const char *karaokefile;
	const char *outfile;
	int windowsize, overlap;
	int delta_min, delta_max;
	int delta_searchmin, delta_searchmax;
	double overdrive;
}
params_t;

void options(params_t *params, int argc, char **argv)
{
	memset(params, 0, sizeof(*params));
	int opt;
	params->windowsize = 4096;
	params->overlap = 512;
	params->delta_min = 0;
	params->delta_max = 0;
	params->delta_searchmin = 4 * 44100;
	params->delta_searchmax = 5 * 44100;
	params->overdrive = 1.2;
	while((opt = getopt(argc, argv, "i:k:o:w:l:d:D:s:S:f:")) != -1)
	{
		switch(opt)
		{
			case 'i':
				params->infile = optarg;
				break;
			case 'k':
				params->karaokefile = optarg;
				break;
			case 'o':
				params->outfile = optarg;
				break;
			case 'l':
				params->overlap = atoi(optarg);
				break;
			case 'w':
				params->windowsize = atoi(optarg);
				break;
			case 'd':
				params->delta_min = atoi(optarg);
				break;
			case 'D':
				params->delta_max = atoi(optarg);
				break;
			case 's':
				params->delta_searchmin = atoi(optarg);
				break;
			case 'S':
				params->delta_searchmax = atoi(optarg);
				break;
			case 'f':
				params->overdrive = atof(optarg);
				break;
			default:
				fprintf(stderr, USAGE, argv[0]);
				exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char **argv)
{
	params_t params;
	options(&params, argc, argv);

	AFfilesetup wavefmt = afNewFileSetup();
	afInitFileFormat(wavefmt, AF_FILE_WAVE);
	AFfilehandle infile = afOpenFile(params.infile, "r", wavefmt);
	if(!infile)
		err(EXIT_FAILURE, "could not open input file");
	AFfilehandle karaokefile = afOpenFile(params.karaokefile, "r", wavefmt);
	if(!karaokefile)
		err(EXIT_FAILURE, "could not open input file");

	{
		// duplicate params
		afInitRate(wavefmt, AF_DEFAULT_TRACK, afGetRate(infile, AF_DEFAULT_TRACK));
		afInitChannels(wavefmt, AF_DEFAULT_TRACK, afGetChannels(infile, AF_DEFAULT_TRACK));
		afInitByteOrder(wavefmt, AF_DEFAULT_TRACK, afGetByteOrder(infile, AF_DEFAULT_TRACK));
		int fmt, wid;
		afGetSampleFormat(infile, AF_DEFAULT_TRACK, &fmt, &wid);
		afInitSampleFormat(wavefmt, AF_DEFAULT_TRACK, fmt, wid);
	}

	AFfilehandle outfile = afOpenFile(params.outfile, "w", wavefmt);

	double rate = afGetRate(infile, AF_DEFAULT_TRACK);
	int channels = afGetChannels(infile, AF_DEFAULT_TRACK);
	afSetVirtualChannels(infile, AF_DEFAULT_TRACK, channels);
	afSetVirtualChannels(karaokefile, AF_DEFAULT_TRACK, channels);
	afSetVirtualChannels(outfile, AF_DEFAULT_TRACK, channels);
	afSetVirtualSampleFormat(infile, AF_DEFAULT_TRACK, AF_SAMPFMT_sample_t, 8 * sizeof(sample_t));
	afSetVirtualSampleFormat(karaokefile, AF_DEFAULT_TRACK, AF_SAMPFMT_sample_t, 8 * sizeof(sample_t));
	afSetVirtualSampleFormat(outfile, AF_DEFAULT_TRACK, AF_SAMPFMT_sample_t, 8 * sizeof(sample_t));

	sample_t *inframes = calloc(channels * (params.windowsize + 2 * params.overlap), sizeof(sample_t));
	sample_t *karaokeframes = calloc(channels * (params.windowsize + 2 * params.overlap), sizeof(sample_t));
	sample_t *outframes = calloc(channels * (params.windowsize + 2 * params.overlap), sizeof(sample_t));

	int nread1, nread2;
	int delta;
	int delta_searchlen = (params.delta_searchmax - params.delta_searchmin) - (params.delta_max - params.delta_min);
	sample_t search_in[channels * (params.delta_searchmax - params.delta_searchmin)];
	sample_t search_kar[channels * (params.delta_searchmax - params.delta_searchmin)];
	afSeekFrame(infile, AF_DEFAULT_TRACK, params.delta_searchmin + params.delta_min);
	afSeekFrame(karaokefile, AF_DEFAULT_TRACK, params.delta_searchmin);
	afReadFrames(infile, AF_DEFAULT_TRACK, search_in, params.delta_searchmax - params.delta_searchmin);
	afReadFrames(karaokefile, AF_DEFAULT_TRACK, search_kar, params.delta_searchmax - params.delta_searchmin);
	double best_score;
	int best_delta;
	int o = params.delta_min;
	for(delta = params.delta_min; delta <= params.delta_max; delta += 64)
	{
		double s = do_fftcompare_samples(delta_searchlen, channels, search_in + channels * (delta - o), search_kar);
		if(delta == params.delta_min || s < best_score)
		{
			best_delta = delta;
			best_score = s;
		}
		fprintf(stderr, "\rSearching (1): %6.2f%%\e[K", (delta - params.delta_min) * 100.0 / (params.delta_max - params.delta_min));
	}
	delta = best_delta;
	fprintf(stderr, "\rIdeal delta: %d samples\e[K\n", best_delta);

	params.delta_min = bound(params.delta_min, delta - 512, params.delta_max);
	params.delta_max = bound(params.delta_min, delta + 512, params.delta_max);
	for(delta = params.delta_min; delta <= params.delta_max; ++delta)
	{
		double s = do_directcompare_samples(delta_searchlen, channels, search_in + channels * (delta - o), search_kar);
		if(delta == params.delta_min || s < best_score)
		{
			best_delta = delta;
			best_score = s;
		}
		fprintf(stderr, "\rSearching (2): %6.2f%%\e[K", (delta - params.delta_min) * 100.0 / (params.delta_max - params.delta_min));
	}
	delta = best_delta;
	fprintf(stderr, "\rIdeal delta: %d samples\e[K\n", delta);

	if(delta > 0)
	{
		afSeekFrame(infile, AF_DEFAULT_TRACK, delta);
		afSeekFrame(karaokefile, AF_DEFAULT_TRACK, 0);
	}
	else
	{
		afSeekFrame(infile, AF_DEFAULT_TRACK, 0);
		afSeekFrame(karaokefile, AF_DEFAULT_TRACK, -delta);
	}


	AFframecount t = afGetFrameCount(infile, AF_DEFAULT_TRACK);
	while((nread1 = afReadFrames(infile, AF_DEFAULT_TRACK, inframes + channels * 2 * params.overlap, params.windowsize)) > 0)
	{
		memset(inframes + channels * (2 * params.overlap + nread1), 0, sizeof(sample_t) * channels * (params.windowsize - nread1));
		nread2 = afReadFrames(karaokefile, AF_DEFAULT_TRACK, karaokeframes + channels * 2 * params.overlap, params.windowsize);
		memset(karaokeframes + channels * (2 * params.overlap + nread2), 0, sizeof(sample_t) * channels * (params.windowsize - nread2));
		memset(outframes + channels * (2 * params.overlap), 0, sizeof(sample_t) * channels * params.windowsize);
		do_filter(params.windowsize, params.overlap, channels, rate, params.overdrive, inframes, karaokeframes, outframes);
		afWriteFrames(outfile, AF_DEFAULT_TRACK, outframes, params.windowsize);
			// TODO this causes an overlap-sized phase shift of output... fix that later by offsetting, or reading the first frame in full
		memmove(inframes, inframes + channels * params.windowsize, sizeof(sample_t) * channels * 2 * params.overlap);
		memmove(karaokeframes, karaokeframes + channels * params.windowsize, sizeof(sample_t) * channels * 2 * params.overlap);
		memmove(outframes, outframes + channels * params.windowsize, sizeof(sample_t) * channels * 2 * params.overlap);
		AFframecount p = afTellFrame(infile, AF_DEFAULT_TRACK);
		fprintf(stderr, "\rProcessing: %6.2f%%\e[K", p * 100.0 / t);
	}
	fprintf(stderr, "\rDone.\e[K\n");
	
	afCloseFile(outfile);
	afCloseFile(karaokefile);
	afCloseFile(infile);

	return 0;
}
