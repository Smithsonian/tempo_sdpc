/** @file daemon.h
 *  @author John C. Houck <jhouck@cfa.harvard.edu>
 *  @date  Oct 2016
 *  @brief Interface for daemon process initialization
 */

#ifndef __DAEMONIZE_INCLUDE__
#define __DAEMONIZE_INCLUDE__ 1

/** Reconfigure the current process as a daemon
 * @param[in] appname   The name of the daemon process 
 * @param[in] logfile_path   Path to a log file for the daemon
 * @return 0 on success, -1 on failure.
 */ 
extern int daemonize (const char *appname, const char *logfile_path);

#endif
