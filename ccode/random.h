/*!
  \file random.h
  \brief Typedefs and defines for the RNG.
*/

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#include <stddef.h> 				  /* size_t */
int ran_init(unsigned int);
int ran_init_(unsigned int seed);
double ran(void);
char *ran_getstate_expert(size_t *len);
char *ran_getstate(void);
int ran_setstate(char *state);
double ran_counter(void);
double ran_(void);
double ran_counter(void);

__END_DECLS
