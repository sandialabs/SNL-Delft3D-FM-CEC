/******************************************************************************
 *
 * Component: OGDI Driver Support Library
 * Purpose: Generic SQL WHERE Expression Evaluator Declarations.
 * Author: Frank Warmerdam <warmerdam@pobox.com>
 *
 ******************************************************************************
 * Copyright (C) 2001 Information Interoperability Institute (3i)
 * Copyright (c) 2010-2013, Even Rouault <even dot rouault at spatialys.com>
 * Permission to use, copy, modify and distribute this software and
 * its documentation for any purpose and without fee is hereby granted,
 * provided that the above copyright notice appear in all copies, that
 * both the copyright notice and this permission notice appear in
 * supporting documentation, and that the name of 3i not be used
 * in advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.  3i makes no
 * representations about the suitability of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 ****************************************************************************/

#ifndef SWQ_H_INCLUDED_
#define SWQ_H_INCLUDED_

#ifndef DOXYGEN_SKIP

#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_core.h"

#include <vector>
#include <set>

#if defined(_WIN32) && !defined(strcasecmp)
#  define strcasecmp stricmp
#endif

// Used for swq_summary.oSetDistinctValues and oVectorDistinctValues
#define SZ_OGR_NULL  "__OGR_NULL__"

typedef enum {
    SWQ_OR,
    SWQ_AND,
    SWQ_NOT,
    SWQ_EQ,
    SWQ_NE,
    SWQ_GE,
    SWQ_LE,
    SWQ_LT,
    SWQ_GT,
    SWQ_LIKE,
    SWQ_ILIKE,
    SWQ_ISNULL,
    SWQ_IN,
    SWQ_BETWEEN,
    SWQ_ADD,
    SWQ_SUBTRACT,
    SWQ_MULTIPLY,
    SWQ_DIVIDE,
    SWQ_MODULUS,
    SWQ_CONCAT,
    SWQ_SUBSTR,
    SWQ_HSTORE_GET_VALUE,
    SWQ_AVG,
    SWQ_MIN,
    SWQ_MAX,
    SWQ_COUNT,
    SWQ_SUM,
    SWQ_CAST,
    SWQ_CUSTOM_FUNC, /* only if parsing done in bAcceptCustomFuncs mode */
    SWQ_ARGUMENT_LIST /* temporary value only set during parsing and replaced by something else at the end */
} swq_op;

typedef enum {
    SWQ_INTEGER,
    SWQ_INTEGER64,
    SWQ_FLOAT,
    SWQ_STRING,
    SWQ_BOOLEAN,  // integer
    SWQ_DATE,     // string
    SWQ_TIME,     // string
    SWQ_TIMESTAMP,// string
    SWQ_GEOMETRY,
    SWQ_NULL,
    SWQ_OTHER,
    SWQ_ERROR
} swq_field_type;

#define SWQ_IS_INTEGER(x) ((x) == SWQ_INTEGER || (x) == SWQ_INTEGER64)

typedef enum {
    SNT_CONSTANT,
    SNT_COLUMN,
    SNT_OPERATION
} swq_node_type;

class swq_field_list;
class swq_expr_node;
class swq_select;
class OGRGeometry;

typedef swq_expr_node *(*swq_field_fetcher)( swq_expr_node *op,
                                             void *record_handle );
typedef swq_expr_node *(*swq_op_evaluator)(swq_expr_node *op,
                                           swq_expr_node **sub_field_values );
typedef swq_field_type (*swq_op_checker)( swq_expr_node *op,
                                          int bAllowMismatchTypeOnFieldComparison );

class swq_custom_func_registrar;

class CPL_UNSTABLE_API swq_expr_node {

    CPL_DISALLOW_COPY_ASSIGN(swq_expr_node)
    swq_expr_node* Evaluate( swq_field_fetcher pfnFetcher,
                             void *record, int nRecLevel );
public:
    swq_expr_node();

    explicit swq_expr_node( const char * );
    explicit swq_expr_node( int );
    explicit swq_expr_node( GIntBig );
    explicit swq_expr_node( double );
    explicit swq_expr_node( OGRGeometry* );
    explicit swq_expr_node( swq_op );

    ~swq_expr_node();

    void           MarkAsTimestamp();
    CPLString      UnparseOperationFromUnparsedSubExpr(char** apszSubExpr);
    char          *Unparse( swq_field_list *, char chColumnQuote );
    void           Dump( FILE *fp, int depth );
    swq_field_type Check( swq_field_list *, int bAllowFieldsInSecondaryTables,
                          int bAllowMismatchTypeOnFieldComparison,
                          swq_custom_func_registrar* poCustomFuncRegistrar,
                          int depth = 0 );
    swq_expr_node* Evaluate( swq_field_fetcher pfnFetcher,
                             void *record );
    swq_expr_node* Clone();

    void           ReplaceBetweenByGEAndLERecurse();

    swq_node_type eNodeType = SNT_CONSTANT;
    swq_field_type field_type = SWQ_INTEGER;

    /* only for SNT_OPERATION */
    void        PushSubExpression( swq_expr_node * );
    void        ReverseSubExpressions();
    int         nOperation = 0;
    int         nSubExprCount = 0;
    swq_expr_node **papoSubExpr = nullptr;

    /* only for SNT_COLUMN */
    int         field_index = 0;
    int         table_index = 0;
    char        *table_name = nullptr;

    /* only for SNT_CONSTANT */
    int         is_null = false;
    GIntBig     int_value = 0;
    double      float_value = 0.0;
    OGRGeometry *geometry_value = nullptr;

    /* shared by SNT_COLUMN, SNT_CONSTANT and also possibly SNT_OPERATION when */
    /* nOperation == SWQ_CUSTOM_FUNC */
    char        *string_value = nullptr; /* column name when SNT_COLUMN */

    static CPLString   QuoteIfNecessary( const CPLString &, char chQuote = '\'' );
    static CPLString   Quote( const CPLString &, char chQuote = '\'' );
};

typedef struct {
    const char*      pszName;
    swq_op           eOperation;
    swq_op_evaluator pfnEvaluator;
    swq_op_checker   pfnChecker;
} swq_operation;

class CPL_UNSTABLE_API swq_op_registrar {
public:
    static const swq_operation *GetOperator( const char * );
    static const swq_operation *GetOperator( swq_op eOperation );
};

class CPL_UNSTABLE_API swq_custom_func_registrar
{
    public:
        virtual ~swq_custom_func_registrar() {}
        virtual const swq_operation *GetOperator( const char * ) = 0;
};

typedef struct {
    char       *data_source;
    char       *table_name;
    char       *table_alias;
} swq_table_def;

class CPL_UNSTABLE_API swq_field_list {
public:
    int count;
    char **names;
    swq_field_type *types;
    int *table_ids;
    int *ids;

    int table_count;
    swq_table_def *table_defs;
};

class CPL_UNSTABLE_API swq_parse_context {
public:
    swq_parse_context() : nStartToken(0), pszInput(nullptr), pszNext(nullptr),
                          pszLastValid(nullptr), bAcceptCustomFuncs(FALSE),
                          poRoot(nullptr), poCurSelect(nullptr) {}

    int        nStartToken;
    const char *pszInput;
    const char *pszNext;
    const char *pszLastValid;
    int        bAcceptCustomFuncs;

    swq_expr_node *poRoot;

    swq_select    *poCurSelect;
};

/* Compile an SQL WHERE clause into an internal form.  The field_list is
** the list of fields in the target 'table', used to render where into
** field numbers instead of names.
*/
int CPL_UNSTABLE_API swqparse( swq_parse_context *context );
int CPL_UNSTABLE_API swqlex( swq_expr_node **ppNode, swq_parse_context *context );
void CPL_UNSTABLE_API swqerror( swq_parse_context *context, const char *msg );

int CPL_UNSTABLE_API swq_identify_field( const char* table_name,
                                         const char *token, swq_field_list *field_list,
                                         swq_field_type *this_type, int *table_id );

CPLErr CPL_UNSTABLE_API swq_expr_compile( const char *where_clause,
                                          int field_count,
                                          char **field_list,
                                          swq_field_type *field_types,
                                          int bCheck,
                                          swq_custom_func_registrar* poCustomFuncRegistrar,
                                          swq_expr_node **expr_root );

CPLErr CPL_UNSTABLE_API swq_expr_compile2( const char *where_clause,
                                           swq_field_list *field_list,
                                           int bCheck,
                                           swq_custom_func_registrar* poCustomFuncRegistrar,
                                           swq_expr_node **expr_root );

/*
** Evaluation related.
*/
int CPL_UNSTABLE_API swq_test_like( const char *input, const char *pattern );

swq_expr_node CPL_UNSTABLE_API *SWQGeneralEvaluator( swq_expr_node *, swq_expr_node **);
swq_field_type CPL_UNSTABLE_API SWQGeneralChecker( swq_expr_node *node, int bAllowMismatchTypeOnFieldComparison );
swq_expr_node CPL_UNSTABLE_API *SWQCastEvaluator( swq_expr_node *, swq_expr_node **);
swq_field_type CPL_UNSTABLE_API SWQCastChecker( swq_expr_node *node, int bAllowMismatchTypeOnFieldComparison );
const char CPL_UNSTABLE_API    *SWQFieldTypeToString( swq_field_type field_type );

/****************************************************************************/

#define SWQP_ALLOW_UNDEFINED_COL_FUNCS 0x01

#define SWQM_SUMMARY_RECORD  1
#define SWQM_RECORDSET       2
#define SWQM_DISTINCT_LIST   3

typedef enum {
    SWQCF_NONE = 0,
    SWQCF_AVG = SWQ_AVG,
    SWQCF_MIN = SWQ_MIN,
    SWQCF_MAX = SWQ_MAX,
    SWQCF_COUNT = SWQ_COUNT,
    SWQCF_SUM = SWQ_SUM,
    SWQCF_CUSTOM
} swq_col_func;

typedef struct {
    swq_col_func col_func;
    char         *table_name;
    char         *field_name;
    char         *field_alias;
    int          table_index;
    int          field_index;
    swq_field_type field_type;
    swq_field_type target_type;
    OGRFieldSubType target_subtype;
    int          field_length;
    int          field_precision;
    int          distinct_flag;
    OGRwkbGeometryType eGeomType;
    int          nSRID;
    swq_expr_node *expr;
} swq_col_def;

class CPL_UNSTABLE_API swq_summary {
public:
    struct Comparator
    {
        bool    bSortAsc;
        swq_field_type eType;

        Comparator() : bSortAsc(true), eType(SWQ_STRING) {}

        bool    operator() (const CPLString&, const CPLString &) const;
    };

    GIntBig     count = 0;

    std::vector<CPLString>          oVectorDistinctValues{};
    std::set<CPLString, Comparator> oSetDistinctValues{};
    double      sum = 0.0;
    double      min = 0.0;
    double      max = 0.0;
    CPLString   osMin{};
    CPLString   osMax{};
};

typedef struct {
    char *table_name;
    char *field_name;
    int   table_index;
    int   field_index;
    int   ascending_flag;
} swq_order_def;

typedef struct {
    int        secondary_table;
    swq_expr_node  *poExpr;
} swq_join_def;

class CPL_UNSTABLE_API swq_select_parse_options
{
public:
    swq_custom_func_registrar* poCustomFuncRegistrar;
    int                        bAllowFieldsInSecondaryTablesInWhere;
    int                        bAddSecondaryTablesGeometryFields;
    int                        bAlwaysPrefixWithTableName;
    int                        bAllowDistinctOnGeometryField;
    int                        bAllowDistinctOnMultipleFields;

                    swq_select_parse_options(): poCustomFuncRegistrar(nullptr),
                                                bAllowFieldsInSecondaryTablesInWhere(FALSE),
                                                bAddSecondaryTablesGeometryFields(FALSE),
                                                bAlwaysPrefixWithTableName(FALSE),
                                                bAllowDistinctOnGeometryField(FALSE),
                                                bAllowDistinctOnMultipleFields(FALSE) {}
};

class CPL_UNSTABLE_API swq_select
{
    void        postpreparse();

    CPL_DISALLOW_COPY_ASSIGN(swq_select)

public:
    swq_select();
    ~swq_select();

    int         query_mode = 0;

    char        *raw_select = nullptr;

    int         PushField( swq_expr_node *poExpr, const char *pszAlias=nullptr,
                           int distinct_flag = FALSE );
    int         result_columns = 0;
    swq_col_def *column_defs = nullptr;
    std::vector<swq_summary> column_summary{};

    int         PushTableDef( const char *pszDataSource,
                              const char *pszTableName,
                              const char *pszAlias );
    int         table_count = 0;
    swq_table_def *table_defs = nullptr;

    void        PushJoin( int iSecondaryTable, swq_expr_node* poExpr );
    int         join_count = 0;
    swq_join_def *join_defs = nullptr;

    swq_expr_node *where_expr = nullptr;

    void        PushOrderBy( const char* pszTableName, const char *pszFieldName, int bAscending );
    int         order_specs = 0;
    swq_order_def *order_defs = nullptr;

    void        SetLimit( GIntBig nLimit );
    GIntBig     limit = -1;

    void        SetOffset( GIntBig nOffset );
    GIntBig     offset = 0;

    swq_select *poOtherSelect = nullptr;
    void        PushUnionAll( swq_select* poOtherSelectIn );

    CPLErr      preparse( const char *select_statement,
                          int bAcceptCustomFuncs = FALSE );
    CPLErr      expand_wildcard( swq_field_list *field_list,
                                 int bAlwaysPrefixWithTableName );
    CPLErr      parse( swq_field_list *field_list,
                       swq_select_parse_options* poParseOptions );

    char       *Unparse();
    void        Dump( FILE * );
};

CPLErr CPL_UNSTABLE_API swq_select_parse( swq_select *select_info,
                                          swq_field_list *field_list,
                                          int parse_flags );

const char CPL_UNSTABLE_API *swq_select_summarize( swq_select *select_info,
                                                   int dest_column,
                                                   const char *value );

int CPL_UNSTABLE_API swq_is_reserved_keyword(const char* pszStr);

char CPL_UNSTABLE_API *OGRHStoreGetValue(const char* pszHStore, 
                                         const char* pszSearchedKey);

#endif /* #ifndef DOXYGEN_SKIP */

#endif /* def SWQ_H_INCLUDED_ */
