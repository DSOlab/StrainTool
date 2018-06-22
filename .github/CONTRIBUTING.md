## Contributing

1. Create an issue and describe your idea
2. [Fork it](https://github.com/demanasta/coulomb2gmt/network#fork-destination-box)
3. Create your feature branch (`git checkout -b my-new-idea`)
4. Commit your changes (`git commit -am 'Add some feature'`)
5. Publish the branch (`git push origin my-new-idea`)
6. Create a new Pull Request
7. Profit! :white_check_mark:

### git structure

* __branch__ master/develop: stucture of main branch
  	* bash scripts
  		* coulomb2gmt.sh: main script
  		* mvclbfiles.sh: assistant script, move and rename files
  		* default-param: configure parameters
  	* functions: bash functions called from main script
  		* messages.sh: help function and ptint messages.
  		* gen_func.sh: general functions.
  		* clbplots.sh: functions for gmt plots
  		* checknum.sh: check number function.
  	* docs: MarkDown templates for issues, pull requests, contributions etc.
  		
  	
* __branch__ documents :
	* tutorial: reference and user guide, tex files.
	* examples: presentation of examples, tex/beamer files.
	
* __branch__ testcase: Include configured files for testing the script.

### Simple guidelines for coding

* __general__
	* Use 80 characters long line.
	* Surround variables with `{}` and use `"${}"` in `if` case.
	* Use comments in coding.
* __gmt functions__
	* Use `-K` `-O` `-V${VRBLEVM}` at the end of each function. 
	* Create a function if a part of the script will be used more than two times.
	* Add printed comments and debug messages in the code.