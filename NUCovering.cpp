#include "fragmentation.h"
#include <locale.h>
#include <chrono>
#include <iostream>
#include <cilk/cilk_api.h>

using namespace std;
using namespace chrono;
/// параметры начальной прямоугольной области
const double g_l1_max = 12.0;
const double g_l2_max = g_l1_max;
const double g_l1_min = 8.0;
const double g_l2_min = g_l1_min;
const double g_l0 = 5.0;

/// точность аппроксимации рабочего пространства
const double g_precision = 0.25;


int main()
{
	setlocale(LC_ALL, "Rus");

	high_resolution_clock::time_point start, end;
	start = high_resolution_clock::now(); //засекли время
	__cilkrts_end_cilk();
	__cilkrts_set_param("nworkers", "4");
	high_level_analysis main_object(g_l0 *(-2), g_l0 * 0, g_l1_min * 3, g_l1_min * 2);

	main_object.GetSolution();
	end = high_resolution_clock::now();
	duration<double> duration = (end - start);
	
	cout << "Количество вычилителей: " << __cilkrts_get_nworkers() << endl;
	cout << "Время исполнения кода: " << duration.count() << " секунд" << std::endl;

	// Внимание! здесь необходимо определить пути до выходных файлов!
	const char* out_files[3] = {	"solution.txt", "nsolution.txt", "boundary.txt" };
	WriteResults( out_files );

	system("pause");
	return 0;
}
