//
// File generated by rootcint at Tue May 21 12:01:48 2019

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME TFileStringDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TFileStringDict.h"

#include "TCollectionProxyInfo.h"
#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void TFileString_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_TFileString(void *p = 0);
   static void *newArray_TFileString(Long_t size, void *p);
   static void delete_TFileString(void *p);
   static void deleteArray_TFileString(void *p);
   static void destruct_TFileString(void *p);
   static void streamer_TFileString(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TFileString*)
   {
      ::TFileString *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TFileString >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TFileString", ::TFileString::Class_Version(), "./TFileString.h", 26,
                  typeid(::TFileString), DefineBehavior(ptr, ptr),
                  &::TFileString::Dictionary, isa_proxy, 0,
                  sizeof(::TFileString) );
      instance.SetNew(&new_TFileString);
      instance.SetNewArray(&newArray_TFileString);
      instance.SetDelete(&delete_TFileString);
      instance.SetDeleteArray(&deleteArray_TFileString);
      instance.SetDestructor(&destruct_TFileString);
      instance.SetStreamerFunc(&streamer_TFileString);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TFileString*)
   {
      return GenerateInitInstanceLocal((::TFileString*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TFileString*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *TFileString::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *TFileString::Class_Name()
{
   return "TFileString";
}

//______________________________________________________________________________
const char *TFileString::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TFileString*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TFileString::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TFileString*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TFileString::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TFileString*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TFileString::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TFileString*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void TFileString::Streamer(TBuffer &R__b)
{
   // Stream an object of class TFileString.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObjString::Streamer(R__b);
      {
         vector<std::string> &R__stl =  mstring;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            string R__t;
            {TString R__str;
             R__str.Streamer(R__b);
             R__t = R__str.Data();}
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, TFileString::IsA());
   } else {
      R__c = R__b.WriteVersion(TFileString::IsA(), kTRUE);
      TObjString::Streamer(R__b);
      {
         vector<std::string> &R__stl =  mstring;
         int R__n=(&R__stl) ? int(R__stl.size()) : 0;
         R__b << R__n;
         if(R__n) {
            vector<std::string>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            {TString R__str((*R__k).c_str());
             R__str.Streamer(R__b);};
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TFileString::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class TFileString.
      TClass *R__cl = ::TFileString::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "mstring", (void*)&mstring);
      R__insp.InspectMember("vector<std::string>", (void*)&mstring, "mstring.", false);
      TObjString::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TFileString(void *p) {
      return  p ? new(p) ::TFileString : new ::TFileString;
   }
   static void *newArray_TFileString(Long_t nElements, void *p) {
      return p ? new(p) ::TFileString[nElements] : new ::TFileString[nElements];
   }
   // Wrapper around operator delete
   static void delete_TFileString(void *p) {
      delete ((::TFileString*)p);
   }
   static void deleteArray_TFileString(void *p) {
      delete [] ((::TFileString*)p);
   }
   static void destruct_TFileString(void *p) {
      typedef ::TFileString current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TFileString(TBuffer &buf, void *obj) {
      ((::TFileString*)obj)->::TFileString::Streamer(buf);
   }
} // end of namespace ROOT for class ::TFileString

namespace ROOT {
   void vectorlEstringgR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void vectorlEstringgR_Dictionary();
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>),0);
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "prec_stl/vector", 49,
                  typeid(vector<string>), DefineBehavior(ptr, ptr),
                  0, &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void vectorlEstringgR_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

/********************************************************
* TFileStringDict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableTFileStringDict();

extern "C" void G__set_cpp_environmentTFileStringDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("TFileString.h");
  G__cpp_reset_tagtableTFileStringDict();
}
#include <new>
extern "C" int G__cpp_dllrevTFileStringDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TFileString */
static int G__TFileStringDict_184_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TFileString* p = NULL;
   char* gvp = (char*) G__getgvp();
   switch (libp->paran) {
   case 1:
     //m: 1
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new TFileString(*(string*) libp->para[0].ref);
     } else {
       p = new((void*) gvp) TFileString(*(string*) libp->para[0].ref);
     }
     break;
   case 0:
     int n = G__getaryconstruct();
     if (n) {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new TFileString[n];
       } else {
         p = new((void*) gvp) TFileString[n];
       }
     } else {
       if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
         p = new TFileString;
       } else {
         p = new((void*) gvp) TFileString;
       }
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TFileStringDictLN_TFileString));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TFileString* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 2
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new TFileString(*(string*) libp->para[0].ref, *(string*) libp->para[1].ref);
   } else {
     p = new((void*) gvp) TFileString(*(string*) libp->para[0].ref, *(string*) libp->para[1].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TFileStringDictLN_TFileString));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         const vector<std::string>& obj = ((TFileString*) G__getstructoffset())->tags();
         result7->ref = (long) (&obj);
         result7->obj.i = (long) (&obj);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         const vector<std::string>& obj = ((const TFileString*) G__getstructoffset())->tags();
         result7->ref = (long) (&obj);
         result7->obj.i = (long) (&obj);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         const string* pobj;
         const string xobj = ((const TFileString*) G__getstructoffset())->name();
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         const string& obj = ((TFileString*) G__getstructoffset())->operator[]((int) G__int(libp->para[0]));
         result7->ref = (long) (&obj);
         result7->obj.i = (long) (&obj);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         const string* pobj;
         const string xobj = ((const TFileString*) G__getstructoffset())->operator[]((int) G__int(libp->para[0]));
         pobj = new string(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 107, (long) ((const TFileString*) G__getstructoffset())->size());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TFileString*) G__getstructoffset())->push_back(*(string*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TFileString*) G__getstructoffset())->add(*(string*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TFileString::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TFileString::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TFileString::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TFileString::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TFileString*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TFileString::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TFileString::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TFileString::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TFileStringDict_184_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TFileString::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__TFileStringDict_184_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   TFileString* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new TFileString(*(TFileString*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__TFileStringDictLN_TFileString));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TFileString G__TTFileString;
static int G__TFileStringDict_184_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (TFileString*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TFileString*) (soff+(sizeof(TFileString)*i)))->~G__TTFileString();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (TFileString*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((TFileString*) (soff))->~G__TTFileString();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__TFileStringDict_184_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TFileString* dest = (TFileString*) G__getstructoffset();
   *dest = *(TFileString*) libp->para[0].ref;
   const TFileString& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TFileString */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTFileStringDict {
 public:
  G__Sizep2memfuncTFileStringDict(): p(&G__Sizep2memfuncTFileStringDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTFileStringDict::*p)();
};

size_t G__get_sizep2memfuncTFileStringDict()
{
  G__Sizep2memfuncTFileStringDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceTFileStringDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TFileStringDictLN_TFileString))) {
     TFileString *G__Lderived;
     G__Lderived=(TFileString*)0x1000;
     {
       TObjString *G__Lpbase=(TObjString*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TFileStringDictLN_TFileString),G__get_linked_tagnum(&G__TFileStringDictLN_TObjString),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TFileStringDictLN_TFileString),G__get_linked_tagnum(&G__TFileStringDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,0);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTFileStringDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TFileStringDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TFileStringDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TFileStringDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TFileStringDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<std::string>",117,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<string>",117,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TFileStringDictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TFileStringDictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<string>",117,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TFileString */
static void G__setup_memvarTFileString(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TFileStringDictLN_TFileString));
   { TFileString *p; p=(TFileString*)0x1000; if (p) { }
   G__memvar_setup((void*)0,117,0,0,G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR),G__defined_typename("vector<std::string>"),-1,4,"mstring=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__TFileStringDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTFileStringDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncTFileString(void) {
   /* TFileString */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TFileStringDictLN_TFileString));
   G__memfunc_setup("TFileString",1099,G__TFileStringDict_184_0_1, 105, G__get_linked_tagnum(&G__TFileStringDictLN_TFileString), -1, 0, 1, 1, 1, 0, "u 'string' - 11 '\"\"' name", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TFileString",1099,G__TFileStringDict_184_0_2, 105, G__get_linked_tagnum(&G__TFileStringDictLN_TFileString), -1, 0, 2, 1, 1, 0, 
"u 'string' - 11 - name u 'string' - 11 - tag", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("tags",431,G__TFileStringDict_184_0_3, 117, G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR), G__defined_typename("vector<std::string>"), 1, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("tags",431,G__TFileStringDict_184_0_4, 117, G__get_linked_tagnum(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR), G__defined_typename("vector<std::string>"), 1, 0, 1, 1, 9, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("name",417,G__TFileStringDict_184_0_5, 117, G__get_linked_tagnum(&G__TFileStringDictLN_string), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("operator[]",1060,G__TFileStringDict_184_0_6, 117, G__get_linked_tagnum(&G__TFileStringDictLN_string), -1, 1, 1, 1, 1, 0, "i - - 0 - i", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("operator[]",1060,G__TFileStringDict_184_0_7, 117, G__get_linked_tagnum(&G__TFileStringDictLN_string), -1, 0, 1, 1, 1, 8, "i - - 0 - i", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("size",443,G__TFileStringDict_184_0_8, 107, -1, G__defined_typename("size_t"), 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("push_back",944,G__TFileStringDict_184_0_9, 121, -1, -1, 0, 1, 1, 1, 0, "u 'string' - 11 - s", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("add",297,G__TFileStringDict_184_0_10, 121, -1, -1, 0, 1, 1, 1, 0, "u 'string' - 11 - s", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TFileStringDict_184_0_11, 85, G__get_linked_tagnum(&G__TFileStringDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&TFileString::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__TFileStringDict_184_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TFileString::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__TFileStringDict_184_0_13, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&TFileString::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__TFileStringDict_184_0_14, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&TFileString::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__TFileStringDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TFileStringDict_184_0_18, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TFileStringDict_184_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TFileString::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TFileStringDict_184_0_20, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TFileString::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__TFileStringDict_184_0_21, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&TFileString::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TFileStringDict_184_0_22, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&TFileString::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("TFileString", 1099, G__TFileStringDict_184_0_23, (int) ('i'), G__get_linked_tagnum(&G__TFileStringDictLN_TFileString), -1, 0, 1, 1, 1, 0, "u 'TFileString' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~TFileString", 1225, G__TFileStringDict_184_0_24, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__TFileStringDict_184_0_25, (int) ('u'), G__get_linked_tagnum(&G__TFileStringDictLN_TFileString), -1, 1, 1, 1, 1, 0, "u 'TFileString' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTFileStringDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalTFileStringDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
}

static void G__cpp_setup_func13() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcTFileStringDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
  G__cpp_setup_func13();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__TFileStringDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_string = { "string" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_TObjString = { "TObjString" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_TFileString = { "TFileString" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR = { "vector<string,allocator<string> >" , 99 , -1 };
G__linked_taginfo G__TFileStringDictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<string,allocator<string> >::iterator>" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTFileStringDict() {
  G__TFileStringDictLN_TClass.tagnum = -1 ;
  G__TFileStringDictLN_TBuffer.tagnum = -1 ;
  G__TFileStringDictLN_TMemberInspector.tagnum = -1 ;
  G__TFileStringDictLN_TObject.tagnum = -1 ;
  G__TFileStringDictLN_string.tagnum = -1 ;
  G__TFileStringDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__TFileStringDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TFileStringDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__TFileStringDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TFileStringDictLN_TObjString.tagnum = -1 ;
  G__TFileStringDictLN_TFileString.tagnum = -1 ;
  G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR.tagnum = -1 ;
  G__TFileStringDictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTFileStringDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_string);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_TObjString);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__TFileStringDictLN_TFileString),sizeof(TFileString),-1,61696,(char*)NULL,G__setup_memvarTFileString,G__setup_memfuncTFileString);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_vectorlEstringcOallocatorlEstringgRsPgR);
   G__get_linked_tagnum_fwd(&G__TFileStringDictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR);
}
extern "C" void G__cpp_setupTFileStringDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTFileStringDict()");
  G__set_cpp_environmentTFileStringDict();
  G__cpp_setup_tagtableTFileStringDict();

  G__cpp_setup_inheritanceTFileStringDict();

  G__cpp_setup_typetableTFileStringDict();

  G__cpp_setup_memvarTFileStringDict();

  G__cpp_setup_memfuncTFileStringDict();
  G__cpp_setup_globalTFileStringDict();
  G__cpp_setup_funcTFileStringDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTFileStringDict();
  return;
}
class G__cpp_setup_initTFileStringDict {
  public:
    G__cpp_setup_initTFileStringDict() { G__add_setup_func("TFileStringDict",(G__incsetup)(&G__cpp_setupTFileStringDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTFileStringDict() { G__remove_setup_func("TFileStringDict"); }
};
G__cpp_setup_initTFileStringDict G__cpp_setup_initializerTFileStringDict;

