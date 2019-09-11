#ifndef MAC_BIN_STORED_HH
#define MAC_BIN_STORED_HH

#include <MAC_TransferExp.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <boolVector.hh>

class MAC_Vector ;

// Tools to perform binary storage.

class MAC_BinStored : public MAC_TransferExp
{
   public: //----------------------------------------------------------

   //-- Access
      
      virtual MAC_Data::Type data_type( void ) const ;

  //-- Binary file related access

      // Is type supported for binary storage through MAC_BinStored ?
      static bool is_type_supported( MAC_Data::Type a_type ) ;

      // Initialize file for writing
      static void init_binary_file( std::string const& file_name ) ;

      // Add `data' to binary file `file_name' and return
      // refering expression.
      static MAC_BinStored const* create_reference(
                    MAC_Object* a_owner,
                    MAC_Data const* a_data,
                    std::string const& file_name,
                    bool local_reference_file_name = false ) ;
      
      // Is `file_name' corresponding to a valid binary file ?
      static bool is_valid_binary_file( std::string const& file_name ) ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      struct BinaryRecord 
      {
            MAC_Data::Type type ;
            size_t length ;
            size_t number ;
            size_t foo[8] ;
      } ;
      
     ~MAC_BinStored( void ) ;
      MAC_BinStored( MAC_BinStored const& other ) ;
      MAC_BinStored& operator=( MAC_BinStored const& other ) ;
      
      MAC_BinStored( MAC_Object* a_owner,
                     MAC_Sequence const* argument_list ) ;
      
      // Retrieve data from binary file `file_name' whose identifier is
      // `record_number'.
      static MAC_Data const* restore_from_binary(
         MAC_Object * a_owner,
         std::string const& file_name,
         size_t record_number ) ;

      // Bad record number
      static const size_t bad_record ;

      // Last record number of file named `file_name' or bad_record if
      // file is an empty valid binary file.
      static size_t last_record_number( std::string const& file_name ) ;
      
   //-- Plug in

      MAC_BinStored( void ) ;
      
      virtual MAC_BinStored* create_replica(
         MAC_Object * a_owner,
         MAC_Sequence const* argument_list ) const ;
      
   //-- Identification

      virtual std::string const& usage( void ) const ;
      
   //-- Arguments

      virtual bool valid_arguments(
                           MAC_Sequence const* some_arguments ) const ;
      
   //-- Transfer implementation
      
      MAC_Data const* data( MAC_Context const* ct ) const ;
      
   //-- Class attributes
      
      static MAC_BinStored const* static_OpComponent ;
      static size_t preamble_length ;
      static size_t last_record ;
      static std::string last_file ;
      
   //-- Attributes
      
      mutable MAC_Data const* my_data ;
};

#endif
