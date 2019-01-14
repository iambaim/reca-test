#!/usr/local/bin/perl

print "Build binary files for eca library\n";

# Delete old files
system "make clean";

# Compile files
system "make caa_main_model1";
system "make caa_main_model2";
system "make caa_main_predict";

# Change 
system "chmod a+x caa_main_model1";
system "chmod a+x caa_main_model2";
system "chmod a+x caa_main_predict";
    
# Find current directory
use Cwd qw();
my $path = Cwd::cwd();

# Copy files to test folder
print "ONLY COPY FILES FOR TESTING!!\n";
print "Current directory: $path\n";
print "copy files to ../../scripts/ \n";
system "cp caa_main_model1 ../../scripts/";
system "cp caa_main_model2 ../../scripts/";
system "cp caa_main_predict ../../scripts/";

# Move files to library folder
system "mv caa_main_model1 ../inst/bin/";
system "mv caa_main_model2 ../inst/bin/";
system "mv caa_main_predict ../inst/bin/";

system "rm *.o";

print "Done!\n";
