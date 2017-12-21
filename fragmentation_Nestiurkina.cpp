#include "fragmentation.h"
#include <fstream>


using namespace std;

cilk::reducer<cilk::op_vector<Box>>solution;
cilk::reducer<cilk::op_vector<Box>>not_solution;
cilk::reducer<cilk::op_vector<Box>>boundary;
cilk::reducer<cilk::op_vector<Box>>temporary_boxes;

using namespace std;
/// функции gj()

//------------------------------------------------------------------------------------------
double g1(double x1, double x2)
{
	return (x1*x1 + x2*x2 - g_l1_max*g_l1_max);
}

//------------------------------------------------------------------------------------------
double g2(double x1, double x2)
{
	return (g_l1_min*g_l1_min - x1*x1 - x2*x2);
}
//------------------------------------------------------------------------------------------
double g3(double x1, double x2)
{
	return (x1*x1 + x2*x2 - g_l2_max*g_l2_max);
}
//------------------------------------------------------------------------------------------
double g4(double x1, double x2)
{
	return (g_l2_min*g_l2_min - x1*x1 - x2*x2);
}
//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(double& min_x, double& min_y, double& x_width, double& y_height)
{
	current_box = Box(min_x, min_y, x_width, y_height);
}
//------------------------------------------------------------------------------------------
low_level_fragmentation::low_level_fragmentation(const Box& box)
{
	current_box = box;
}
//------------------------------------------------------------------------------------------
void low_level_fragmentation::VerticalSplitter(const Box& box, boxes_pair& vertical_splitter_pair)
{
	// Определили функцию
	double x, y, width, height;
	box.GetParameters(x, y, width, height);
	vertical_splitter_pair.first = Box(x, y, width / 2, height);
	vertical_splitter_pair.second = Box(x + width / 2, y, width / 2, height);	
}
//------------------------------------------------------------------------------------------
void low_level_fragmentation::HorizontalSplitter(const Box& box, boxes_pair& horizontal_splitter_pair)
{
	// Определили функцию
	double x, y, width, height;
	box.GetParameters(x, y, width, height);
	horizontal_splitter_pair.first = Box(x, y, width, height / 2);
	horizontal_splitter_pair.second = Box(x, y + height / 2, width, height / 2);
}
//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetNewBoxes(const Box& box, boxes_pair& new_pair_of_boxes)
{
	// Определили функцию
	double x, y, width, height;
	box.GetParameters(x, y, width, height);
	if (width >= height)
		VerticalSplitter(box, new_pair_of_boxes);
	else
		HorizontalSplitter(box, new_pair_of_boxes);
}
//------------------------------------------------------------------------------------------
unsigned int low_level_fragmentation::FindTreeDepth()
{
	double box_diagonal = current_box.GetDiagonal();
	if (box_diagonal <= g_precision)
	{
		return 0;
	}
	else
	{
		boxes_pair new_boxes;
		// допустим, разобьем начальную область по ширине
		VerticalSplitter(current_box, new_boxes);
 		unsigned int tree_depth = 1;
		box_diagonal = new_boxes.first.GetDiagonal();

		if (box_diagonal <= g_precision)
		{
			return tree_depth;
		}
		else
		{
			for (;;)
			{
				GetNewBoxes(new_boxes.first, new_boxes);
				++tree_depth;
				box_diagonal = new_boxes.first.GetDiagonal();

				if (box_diagonal <= g_precision)
				{
					break;
				}
			}
			return tree_depth;
		}
	}
}
//------------------------------------------------------------------------------------------
int low_level_fragmentation::ClasifyBox(const min_max_vectors& vects)
{
	// Определили функцию
	int count = 0;
	int firstL = vects.first.size();
	int secondL = vects.second.size();

	if (vects.first[0] == 0 && vects.second[0] == 0) return 3; //Граница прямоугольника

	for (int i = 0; i < 4; i++)
	{
		if (vects.first[i] < 0) count++; //Проверка выполнения условия 4
		if (vects.second[i] > 0) return 0; //Проверка выполнения условия 5
	}

	if (count == 4) return 1;
	else return 2; //Разбиваем на 2 новых		
}
//------------------------------------------------------------------------------------------
void low_level_fragmentation::GetBoxType(const Box& box)
{
	// Определили функцию

	int type;
	min_max_vectors vecs;
	GetMinMax(box, vecs);
	boxes_pair pair;

	type = ClasifyBox(vecs);

	if (type == 0) //Убираем прямоугольник из рассмотрения
		not_solution->push_back(box);

	if (type == 1) //Прямоугольник, в котором есть решение (нужный экстремум)
		solution->push_back(box);

	if (type == 2) { //Прямоугольник, который надо разбить
		GetNewBoxes(box, pair);
		temporary_boxes->push_back(pair.first);
		temporary_boxes->push_back(pair.second);
	}

	if (type == 3) //Граница
		boundary->push_back(box);
}
//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis(double& min_x, double& min_y, double& x_width, double& y_height) :
low_level_fragmentation(min_x, min_y, x_width, y_height) {}
//------------------------------------------------------------------------------------------
high_level_analysis::high_level_analysis(Box& box) : low_level_fragmentation(box) {}
//------------------------------------------------------------------------------------------
void high_level_analysis::GetMinMax(const Box& box, min_max_vectors& min_max_vecs)

{

	std::vector<double> g_min;
	std::vector<double> g_max;

	double a1min, a2min, a1max, a2max;
	double xmin, xmax, ymin, ymax;

	box.GetParameters(xmin, ymin, xmax, ymax);

	xmax = xmin + xmax;
	ymax = ymin + ymax;

	double curr_box_diagonal = box.GetDiagonal();

	if (curr_box_diagonal <= g_precision)
	{
		g_min.push_back(0);
		g_max.push_back(0);

		min_max_vecs.first = g_max;
		min_max_vecs.second = g_min;
		return;
	}

	// MIN
	// функция g1(x1,x2)
	a1min = __min(abs(xmin), abs(xmax));
	a2min = __min(abs(ymin), abs(ymax));
	g_min.push_back(g1(a1min, a2min));

	// функция g2(x1,x2)
	a1min = __max(abs(xmin), abs(xmax));
	a2min = __max(abs(ymin), abs(ymax));
	g_min.push_back(g2(a1min, a2min));

	// функция g3(x1,x2)
	a1min = __min(abs(xmin - g_l0), abs(xmax - g_l0));
	a2min = __min(abs(ymin), abs(ymax));
	g_min.push_back(g3(a1min, a2min));

	// функция g4(x1,x2)
	a1min = __max(abs(xmin - g_l0), abs(xmax - g_l0));
	a2min = __max(abs(ymin), abs(ymax));
	g_min.push_back(g4(a1min, a2min));

	// MAX
	// функция g1(x1,x2)
	a1max = __max(abs(xmin), abs(xmax));
	a2max = __max(abs(ymin), abs(ymax));
	g_max.push_back(g1(a1max, a2max));

	// функция g2(x1,x2)
	a1max = __min(abs(xmin), abs(xmax));
	a2max = __min(abs(ymin), abs(ymax));
	g_max.push_back(g2(a1max, a2max));

	// функция g3(x1,x2)
	a1max = __max(abs(xmin - g_l0), abs(xmax - g_l0));
	a2max = __max(abs(ymin), abs(ymax));
	g_max.push_back(g3(a1max, a2max));

	// функция g4(x1,x2)
	a1max = __min(abs(xmin - g_l0), abs(xmax - g_l0));
	a2max = __min(abs(ymin), abs(ymax));
	g_max.push_back(g4(a1max, a2max));

	min_max_vecs.first = g_max;
	min_max_vecs.second = g_min;
}
//------------------------------------------------------------------------------------------
void high_level_analysis::GetSolution()
{
	// Определили функцию
	int length = FindTreeDepth()+1;	
	boxes_pair pair;
	temporary_boxes->push_back(current_box);

	for (int i = 0; i < length; i++) 
	{
		vector<Box> tmp;
		temporary_boxes.move_out(tmp);
		int number_of_box_on_level = tmp.size();
		vector<Box> curr_boxes(tmp);
		tmp.clear();
		temporary_boxes.set_value(tmp);
		cilk_for (int j = 0; j < number_of_box_on_level; j++)
		{ 
			GetBoxType(curr_boxes[j]);
		}
	}
}
//------------------------------------------------------------------------------------------
void WriteResults(const char* file_names[])
{
	// Определии функцию
	double _xmin, _xmax, _w, _h;
	ofstream sol(file_names[0]);
	vector<Box> sltn;
	solution.move_out(sltn);
	sol << "x_min"  << "\t" << "y_min" << "\t" << "width" << "\t" << "height" << endl;
	for (int i = 0; i < solution.size(); i++) {
		sltn[i].GetParameters(_xmin, _xmax, _w, _h);
		sol << _xmin << "\t" << _xmax << "\t" << _w << "\t" << _h << endl;
	}
	sol.close();

	ofstream nsol(file_names[1]);
	vector<Box> nsltn;
	not_solution.move_out(nsltn);
	nsol << "x_min" << "\t" << "y_min" << "\t" << "width" << "\t" << "height" << endl;
	for (int i = 0; i < not_solution.size(); i++) {
		nsltn[i].GetParameters(_xmin, _xmax, _w, _h);
		nsol << _xmin << "\t" << _xmax << "\t" << _w << "\t" << _h << endl;
	}
	nsol.close();

	ofstream bndr(file_names[2]);
	vector<Box> bndr;
	boundary.move_out(bndr);
	bndr << "x_min" << "\t" << "y_min" << "\t" << "width" << "\t" << "height" << endl;
	for (int i = 0; i < boundary.size(); i++) {
		bndr[i].GetParameters(_xmin, _xmax, _w, _h);
		bndr << _xmin << "\t" << _xmax << "\t" << _w << "\t" << _h << endl;
	}
	bndr.close();
}
