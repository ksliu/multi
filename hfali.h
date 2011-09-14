#ifndef  HFALI_H
#define  HFALI_H
#include "protein.h"
#include "corres.h"

struct SFP
{
	int ia, ib, score;
	SFP(int a = 0, int b = 0, int c = 0) :
		ia(a), ib(b), score(c)
	{
	}
};
inline bool deSFPcmp(const SFP & a, const SFP & b)
{
	return a.score > b.score;
}

class hfali
{
public:
	hfali(const char *subject, const char *query);
	~hfali();

	enum IM
	{
		IM_once, IM_zoom_in
	};

	void solve(int *a1, int *a2, int *b1, int *b2, IM method = IM_once);

	void output(std::ostream & fs = std::cout) const;
	void output_script(const char *fn) const;

private:
	void fill(double fract, double abs_err); // fill short-SFP list (1/2, 5/6, 1), aka zoom in
	void tune(double core, double elong);

	void add(int fromBegin, int fromEnd, int toBegin, int toEnd, int *outAlign);
	void move(double rot[][3], double trans[3]) const;

	static void adjacent_window_score(std::string a, std::string b, int width,
			int(*score)(char, char), int lower_socre, std::vector<SFP> &q);

	static void window_score(std::string a, std::string b, int width,
			int(*score)(char, char), int lower_socre, std::vector<SFP> &q);

	static void shave_adjacent_window_score(std::vector<SFP> &q, int d);

	IM method_log;

	protein *a, *b;
	corres cor;

	static const int SW = 8, SC = 0, N=3000; // width of fragment, cutoff of short SFP CLESUM score
	std::vector<SFP> sq; // short (fragment) SFP lists
	std::vector<int> scon; // consistence-index of s-SFP
	int init_ali[N], core_ali[N], ali[N];
};

#endif  /*HFALI_H*/
