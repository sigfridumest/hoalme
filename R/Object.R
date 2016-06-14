# Defines the class Object and is it constructor.
 Object <- function() {
   # Create a new environment and wrap it up as an attribute to the
   # smallest R object available, that is, NA.
   this <- NA;
   attr(this, ".env") <- new.env();
   class(this) <- "Object";
   this;
 }

 # Returns a string representation of the Object. By default it
 # is "{class name}: 0x{memory address}".
 as.character.Object <- function(this) {
   # getAddress() is a private method for extracting the pointer address
   # where the environment situated all to be able to distinguish one
   # Object from another. Note really necessary though.
   getAddress <- function(this) {
     con <- textConnection("pointer", open="w");
     on.exit(close(con));
     sink(con);
     print.default(attr(this, ".env"));
     sink();
     pointer <- substr(pointer[1], 15, 21);
     pointer;
   }

   paste(getClass(this), ": 0x", getAddress(this), sep="");
 }


 # Clone the Object by copying all of its content.
 clone.Object <- function(this) {
   # Copy the reference.
   clone <- this;

   # Create a new environment, i.e. a new Object.
   clone.env <- new.env();
   attr(clone, ".env") <- clone.env;

   this.env <- attr(this, ".env");

   # Copy all variables in the environment.
   for (name in ls(envir=this.env, all.names=TRUE)) {
     value <- get(name, envir=this.env, inherits=FALSE);
     assign(name, value, envir=clone.env);
   }

   clone;
 }

 clone <- function(...) UseMethod("clone");      # New generic function.

 # Get the class name of the Object (not its superclasses).
 getClass.Object <- function(this) {
   class(this)[1];
 }

 getClass <- function(...) UseMethod("getClass");   # New generic function.

 # Print information about the Object to standard output.
 print.Object <- function(this, ...) {
   print(as.character(this));
 }

 # Map ${name} to return the value of the field {name}, which lives
 # inside the private environment.
 "$.Object" <- function(this, name) {
   get(name, envir=attr(this, ".env"), inherits= F);
 }

 "[[.Object" <- function(this, name) {
   get(name, envir=attr(this, ".env"), inherits= F);
 }


 # Map ${name} <- {value} to assign {value} to field {name}, which
 # lives inside the private environment.
 "$<-.Object" <- function(this, name, value) {
   assign(name, value, envir=attr(this, ".env"));
   this;
 }

# Create instance of class {classname}, by taking another Object,
# add {classname} to the class list and add all the named values
# in ... as fields to the new Object.
extend.Object <- function(this, classname, ...) {
  fields <- list(...);
  names <- names(fields);
  for (name in names)
    assign(name, fields[[name]], envir=attr(this, ".env"));
  class(this) <- c(classname, class(this));
  this;
}

extend <- function(...) UseMethod("extend");
