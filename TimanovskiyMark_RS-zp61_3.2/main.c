#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define HARD_ITER_LIMIT 1000000 

double f(double x);
double df(double x);
void clear_stdin(void);
int check_derivative_monotony(double a, double b);
int ask_iteration_limit_action(unsigned long iter, double x_curr, double fx_curr);

double solve_bisection(double a, double b, double eps, unsigned long max_iter,
                       unsigned long *out_iters, int *out_status);
double solve_chord(double a, double b, double eps, unsigned long max_iter,
                   unsigned long *out_iters, int *out_status);
double solve_newton(double a, double b, double x0, double eps, unsigned long max_iter,
                    unsigned long *out_iters, int *out_status);

int main(void)
{
    int method = 0;                 
    double a = 0.0, b = 0.0;       
    double x0 = 0.0;               
    double eps = 1e-6;              
    unsigned long max_iter = 100;   
    unsigned long iters = 0;        
    int status = 0;                
    double root = 0.0;              
    clock_t t_start, t_end;        
    double elapsed_sec = 0.0;       

    printf("============================================================\n");
    printf("  Laboratory Work #3.2 | Solving Nonlinear Equations\n");
    printf("  Variant #5: f(x) = x^3 - 3x^2 - 6x + 8 = 0\n");
    printf("============================================================\n\n");


    do
    {
        printf("Select numerical method:\n");
        printf("  1 - Bisection method (metod polovynnoho dilennia)\n");
        printf("  2 - Chord method (metod khord)\n");
        printf("  3 - Newton's tangent method (metod dotychnykh)\n");
        printf("Your choice (1-3): ");

        if (scanf("%d", &method) != 1 || method < 1 || method > 3)
        {
            printf("Input error! Please enter 1, 2, or 3.\n\n");
            clear_stdin();
            method = 0;
        }
    } while (method == 0);
    clear_stdin();

    while (1)
    {
        printf("\nEnter isolation interval boundaries [a, b] (example: 3.2 4.5): ");
        if (scanf("%lf %lf", &a, &b) != 2)
        {
            printf("Input error! Enter two valid real numbers.\n");
            clear_stdin();
            continue;
        }
        clear_stdin();

        if (a >= b)
        {
            printf("Error: Left boundary 'a' must be strictly less than 'b'.\n");
            continue;
        }


        if (f(a) * f(b) > 0.0)
        {
            printf("Warning: f(a) = %.6lf and f(b) = %.6lf have the same sign!\n", f(a), f(b));
            printf("The interval [%.4lf, %.4lf] may not isolate a single root. Try again.\n", a, b);
            continue;
        }
        break;
    }


    printf("\nChecking derivative f'(x) on [%.4lf, %.4lf]:\n", a, b);
    printf("  f'(a) = %.6lf, f'(b) = %.6lf\n", df(a), df(b));
    if (check_derivative_monotony(a, b))
    {
        printf("  Result: f'(x) preserves its sign -> function is monotonic (unique root confirmed).\n");
    }
    else
    {
        printf("  Warning: f'(x) changes sign on [a, b] -> derivative may vanish inside the interval!\n");
    }


    while (1)
    {
        printf("\nEnter required accuracy eps (e.g., 0.000001): ");
        if (scanf("%lf", &eps) != 1 || eps <= 0.0)
        {
            printf("Input error! Accuracy eps must be positive.\n");
            clear_stdin();
            continue;
        }
        clear_stdin();
        break;
    }

    while (1)
    {
        printf("Enter iteration step limit N_max before prompt (e.g., 5 or 100): ");
        if (scanf("%lu", &max_iter) != 1 || max_iter == 0)
        {
            printf("Input error! N_max must be a positive integer.\n");
            clear_stdin();
            continue;
        }
        clear_stdin();
        break;
    }


    if (method == 3)
    {
        while (1)
        {
            printf("Enter initial guess x0 in [%.4lf, %.4lf] (recommended f(x0)*f''(x0) > 0): ", a, b);
            if (scanf("%lf", &x0) != 1 || x0 < a || x0 > b)
            {
                printf("Input error! x0 must lie within [%.4lf, %.4lf].\n", a, b);
                clear_stdin();
                continue;
            }
            clear_stdin();
            break;
        }
    }


    t_start = clock();

    switch (method)
    {
        case 1:
            root = solve_bisection(a, b, eps, max_iter, &iters, &status);
            break;
        case 2:
            root = solve_chord(a, b, eps, max_iter, &iters, &status);
            break;
        case 3:
            root = solve_newton(a, b, x0, eps, max_iter, &iters, &status);
            break;
    }

    t_end = clock();
    elapsed_sec = (double)(t_end - t_start) / CLOCKS_PER_SEC;


    printf("\n------------------------------------------------------------\n");
    if (status == 1)
    {
        printf("Status          : Root found successfully!\n");
    }
    else if (status == 0)
    {
        printf("Status          : Stopped by user (intermediate result).\n");
    }
    else
    {
        printf("Status          : Calculation error (division by zero or divergence).\n");
    }

    printf("Root value (x)  : %.10lf\n", root);
    printf("Function f(x)   : %.10e\n", f(root));
    printf("Iterations      : %lu\n", iters);
    printf("CPU time        : %.6f sec (%ld clock ticks)\n",
           elapsed_sec, (long)(t_end - t_start));
    printf("------------------------------------------------------------\n");

    return 0;
}


double f(double x)
{
    return x * x * x - 3.0 * x * x - 6.0 * x + 8.0;
}


double df(double x)
{
    return 3.0 * x * x - 6.0 * x - 6.0;
}


void clear_stdin(void)
{
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF)
    {
    }
}


int check_derivative_monotony(double a, double b)
{
    const int steps = 100;
    double h = (b - a) / steps;
    double d_first = df(a);

    if (fabs(d_first) < 1e-12)
    {
        return 0;
    }

    for (int i = 1; i <= steps; i++)
    {
        double d_curr = df(a + i * h);
        if (d_first * d_curr <= 0.0)
        {
            return 0;
        }
    }
    return 1;
}


int ask_iteration_limit_action(unsigned long iter, double x_curr, double fx_curr)
{
    int choice = 0;
    printf("\n[!] Reached %lu iterations. Current approximation: x = %.8lf, f(x) = %.4e\n",
           iter, x_curr, fx_curr);
    printf("Choose next action:\n");
    printf("  1 (a) - Continue with the same number of iterations\n");
    printf("  2 (b) - Run until the root is found\n");
    printf("  3 (c) - Exit program and display intermediate result\n");

    do
    {
        printf("Your choice (1-3): ");
        if (scanf("%d", &choice) != 1 || choice < 1 || choice > 3)
        {
            printf("Invalid choice! Enter 1, 2, or 3.\n");
            clear_stdin();
            choice = 0;
        }
    } while (choice == 0);
    clear_stdin();

    return choice;
}


double solve_bisection(double a, double b, double eps, unsigned long max_iter,
                       unsigned long *out_iters, int *out_status)
{
    double fa = f(a);
    double fb = f(b);
    double c = 0.5 * (a + b);
    double fc = f(c);
    unsigned long iter = 0;
    int run_to_end = 0;

    if (fabs(fa) <= eps)
    {
        *out_iters = 0;
        *out_status = 1;
        return a;
    }
    if (fabs(fb) <= eps)
    {
        *out_iters = 0;
        *out_status = 1;
        return b;
    }

    while (fabs(b - a) / 2.0 > eps && iter < HARD_ITER_LIMIT)
    {
        iter++;
        c = 0.5 * (a + b);
        fc = f(c);

        if (fabs(fc) <= eps)
        {
            *out_iters = iter;
            *out_status = 1;
            return c;
        }

        if (fa * fc < 0.0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }

        if (!run_to_end && (iter % max_iter == 0) && (fabs(b - a) / 2.0 > eps))
        {
            int action = ask_iteration_limit_action(iter, 0.5 * (a + b), f(0.5 * (a + b)));
            if (action == 2)
            {
                run_to_end = 1;
            }
            else if (action == 3)
            {
                *out_iters = iter;
                *out_status = 0;
                return 0.5 * (a + b);
            }
        }
    }

    c = 0.5 * (a + b);
    *out_iters = iter;
    *out_status = 1;
    return c;
}


double solve_chord(double a, double b, double eps, unsigned long max_iter,
                   unsigned long *out_iters, int *out_status)
{
    double fa = f(a);
    double fb = f(b);
    double x = a;
    double x_prev = a;
    double fx = fa;
    unsigned long iter = 0;
    int run_to_end = 0;

    while (iter < HARD_ITER_LIMIT)
    {
        if (fabs(fb - fa) < 1e-15)
        {
            *out_iters = iter;
            *out_status = -1;
            return x;
        }

        iter++;
        x = a - fa * (b - a) / (fb - fa);
        fx = f(x);

        if (fabs(fx) <= eps || (iter > 1 && fabs(x - x_prev) <= eps))
        {
            *out_iters = iter;
            *out_status = 1;
            return x;
        }

        if (fa * fx < 0.0)
        {
            b = x;
            fb = fx;
        }
        else
        {
            a = x;
            fa = fx;
        }

        x_prev = x;

        if (!run_to_end && (iter % max_iter == 0))
        {
            int action = ask_iteration_limit_action(iter, x, fx);
            if (action == 2)
            {
                run_to_end = 1;
            }
            else if (action == 3)
            {
                *out_iters = iter;
                *out_status = 0;
                return x;
            }
        }
    }

    *out_iters = iter;
    *out_status = -1;
    return x;
}


double solve_newton(double a, double b, double x0, double eps, unsigned long max_iter,
                    unsigned long *out_iters, int *out_status)
{
    double x = x0;
    double x_next = x0;
    unsigned long iter = 0;
    int run_to_end = 0;

    (void)a; 
    (void)b;

    while (iter < HARD_ITER_LIMIT)
    {
        double fx = f(x);
        double dfx = df(x);

        if (fabs(dfx) < 1e-15)
        {
            *out_iters = iter;
            *out_status = -1;
            return x;
        }

        iter++;
        x_next = x - fx / dfx;

        if (fabs(x_next - x) <= eps && fabs(f(x_next)) <= eps)
        {
            *out_iters = iter;
            *out_status = 1;
            return x_next;
        }

        x = x_next;

        if (!run_to_end && (iter % max_iter == 0))
        {
            int action = ask_iteration_limit_action(iter, x, f(x));
            if (action == 2)
            {
                run_to_end = 1;
            }
            else if (action == 3)
            {
                *out_iters = iter;
                *out_status = 0;
                return x;
            }
        }
    }

    *out_iters = iter;
    *out_status = -1;
    return x;
}