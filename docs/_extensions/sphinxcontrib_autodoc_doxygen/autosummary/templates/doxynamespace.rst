{% set title = name + ' module reference' %}
{{ '=' * title|length }}
{{ title }}
{{ '=' * title|length }}

.. f:module:: {{ name }}

{% for line in brief_desc %}
{{ line }}
{% endfor %}

.. _DETA{{ name }}:

-----------------------
Detailed Description
-----------------------

{% for line in detailed_desc %}
{{ line }}
{% endfor %}

{% if types %}
-----------------------
Type Documentation
-----------------------

{% for type in types %}
.. f:type:: {{ type.name }}

{% for line in type.brief %}
   {{ line }}
{% endfor %}

{% for field in type.fields %}
   {{ field }}
{% endfor %}

{% endfor %}
{% endif %}
{% if methods %}
-------------------------------------
Function/Subroutine Documentation
-------------------------------------

{% for method in methods %}
.. f:{{ method.directive }}:: {{ method.signame }}{{ method.argsstring }}

{% for line in method.brief %}
   {{ line }}
{% endfor %}

{% for line in method.detailed %}
   {{ line }}
{% endfor %}

{% endfor %}
{% endif %}
