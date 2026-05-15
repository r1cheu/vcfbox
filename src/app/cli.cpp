#include "cli.h"

#include <exception>
#include <iostream>

#include <CLI/CLI.hpp>

#include "app/subcommands.h"

namespace
{
template <typename Function>
int run_command(Function run)
{
    try
    {
        run();
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
}  // namespace

namespace vcfbox
{
int run_cli(int argc, char** argv)
{
    CLI::App app{"A simple cli tool contains some operation on vcf"};
    argv = app.ensure_utf8(argv);
    app.require_subcommand(1);

    CombineOptions combine_options;
    ParentageMatrixOptions parentage_matrix_options;
    ParentageTestOptions parentage_test_options;

    auto* combine_command = add_combine_command(app, combine_options);
    auto parentage = add_parentage_command(
        app, parentage_matrix_options, parentage_test_options);

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    if (*combine_command)
    {
        return run_command([&] { run_combine_command(combine_options); });
    }
    if (*parentage.build_matrix)
    {
        return run_command(
            [&] { build_parentage_matrices(parentage_matrix_options); });
    }
    if (*parentage.test)
    {
        return run_command(
            [&] { test_parentage(parentage_test_options); });
    }

    return 0;
}
}  // namespace vcfbox
